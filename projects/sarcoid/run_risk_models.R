
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "sarcoid"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Connect to db
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

### Load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

num_cores <- 5
cp_selected <- delay_params$cp

#update demo1 and demo2
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
demo2 <- demo2 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")

rural_ids <- rural_visits %>% inner_join(index_cases %>% distinct(patient_id)) %>% distinct(patient_id)


# update all_dx_visits
load(paste0(delay_base_path,"all_dx_visits.RData"))
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

# update sim_obs
sim_obs <- sim_obs %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<0)

# update time map
load(paste0(delay_base_path,"delay_tm.RData"))

problem_patient_ids <- tm %>% filter(days_since_index == 0) %>% distinct(patient_id) %>%
  anti_join(tm %>% distinct(patient_id),.)

tm <- tm %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

#Load in the setting types
load(paste0(delay_base_path,"visit_info.RData"))

# update tm_stdplac
tm_stdprov <- tm_stdprov %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

load(paste0(sim_out_path,"/sim_res_ssd.RData"))
sim_res_ssd <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  select(trial, boot_trial, res) %>% 
  unnest(res) %>% 
  select(-boot_id)

load(paste0(sim_out_path,"sim_obs_reduced.RData"))
sim_obs_reduced <- sim_obs_reduced %>% 
  mutate(days_since_index = -period) %>% 
  select(-period)

sim_res_sim_obs <- sim_res_ssd %>% 
  inner_join(sim_obs_reduced, by = "obs") 

load(paste0(sim_out_path,"boot_data.RData"))
boot_data <- boot_data %>% select(boot_trial, boot_sample)

n_trials <- nrow(distinct(sim_res_ssd,trial))

# Demo data
reg_demo <- demo1 %>% 
  mutate(female=(sex==2),
         age = index_year-dobyr,
         stdrace = as.numeric(stdrace)) %>% 
  left_join(tibble(stdrace = c(0,1,2,4,9),
                   race = c("Missing/Unknown","White","Black","Hispanic","Other")),
            by = "stdrace") %>% 
  mutate(race = fct_relevel(race,"White"))

reg_demo <- reg_demo %>% 
  left_join(demo2 %>% 
              filter(index_date<=dtend & index_date>=dtstart) %>% 
              mutate(msa_new = msa %in% c("0","")) %>% 
              mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
              mutate(source = as.factor(source)) %>% 
              distinct(patient_id,source,msa=msa_new),
            by = "patient_id")

age_cats <- c(-1,17,34,44,54,64,130)

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,age,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))

# load drug codes

load("/Shared/Statepi_Diagnosis/projects/sarcoid/code_sets/drug_risk_factor_codes.RData")

rx_visits <- con %>% tbl("all_rx_visits") %>% 
  filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
  collect()

rx_visits <- rx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))

gc()

## Prepare Oral Steroids indicators -----------------------------------------------

oral_steroids_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% unique(steroids$ndcnum))

reg_demo <- reg_demo %>% 
  left_join(oral_steroids_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(oral_steroids_window=1L)) %>% 
  mutate_at(vars(oral_steroids_window),~replace_na(.,0L))

## Prepare Inhalers indicators -----------------------------------------------

inhalers_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% unique(inhalers$ndcnum))

reg_demo <- reg_demo %>% 
  left_join(inhalers_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(inhalers_window=1L)) %>% 
  mutate_at(vars(inhalers_window),~replace_na(.,0L))

## Prepare Other Drug indicators -----------------------------------------------
other_drugs <- other_drugs %>% mutate(drug_class = stringr::str_replace_all(drug_class, " ", "_" ))
other_drugs <- other_drugs %>% 
  filter(!drug_class %in% c("antifungal", "antiinflammatory")) %>% 
  mutate(drug_class = ifelse(drug_class %in% c("antiacid", "ppi"), "antiacid_ppi", drug_class))

other_drugs_rx_vis <- rx_visits %>% 
  inner_join(other_drugs) %>% 
  select(-ndcnum, -daysupp) %>% 
  distinct() %>% 
  mutate(ind = 1) %>% 
  pivot_wider(id_cols = c(patient_id, date, index_date),
              values_from = ind,
              names_from = drug_class) %>% 
  mutate_at(vars(antiacid_ppi:cough_suppressant),~replace_na(.,0L))

reg_demo <- reg_demo %>% 
  left_join(other_drugs_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              rename_with(~str_c(.,"_drugs_window"), antiacid_ppi:cough_suppressant) %>% 
              group_by(patient_id) %>% 
              summarise(across(antiacid_ppi_drugs_window:cough_suppressant_drugs_window, max)) %>% 
              ungroup()) %>% 
  mutate_at(vars(antiacid_ppi_drugs_window:cough_suppressant_drugs_window),~replace_na(.,0L))

## Add oral sterioids, inhalers, and other drugs to time map
oral_steroids_rx_vis <- oral_steroids_rx_vis %>% distinct(patient_id, admdate = date)
inhalers_rx_vis <- inhalers_rx_vis %>% distinct(patient_id, admdate = date)
other_drugs_rx_vis <- other_drugs_rx_vis %>% select(-index_date) %>% rename(admdate = date)

tm <- tm %>% rename(admdate = svcdate) %>% 
  left_join(oral_steroids_rx_vis %>% mutate(oral_steroids = 1L), by = c("patient_id", "admdate")) %>% 
  left_join(inhalers_rx_vis %>% mutate(inhalers = 1L), by = c("patient_id", "admdate")) %>% 
  left_join(other_drugs_rx_vis, by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(oral_steroids, inhalers, antiacid_ppi:cough_suppressant),~replace_na(.,0L)) 


## Add morbid obesity and any obesity to reg_demo
morbid_obesity_codes <- tibble(dx = c("27801", "E6601", "E662"),
                               dx_ver = c(9, rep(10, 2)))
any_obesity <-  tibble(dx =c("27800", "27801", "27802", "27803", "E66", "E660", "E6601", "E6609", "E661", "E662", "E663", "E668", "E669"),
                       dx_ver = c(rep(9, 4), rep(10, 9)))

con %>% 
  copy_to(any_obesity, temporary = T) 

obesity_dx_visits <- con %>% 
  tbl("all_dx_visits") %>%
  filter(days_since_index<0) %>% 
  inner_join(tbl(con, "any_obesity")) %>% 
  collect()

obesity_inds <- obesity_dx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  mutate(any_obesity = 1) %>% 
  left_join(morbid_obesity_codes %>% mutate(morbid_obesity = 1)) %>% 
  mutate(morbid_obesity = replace_na(morbid_obesity,0)) %>% 
  group_by(patient_id) %>% 
  summarise(any_obesity = max(any_obesity),
            morbid_obesity = max(morbid_obesity)) %>% 
  ungroup() %>% 
  select(-morbid_obesity)

reg_demo <- reg_demo %>%
  left_join(obesity_inds) %>% 
  mutate_at(vars(any_obesity),~replace_na(.,0L))

############################
#### Prepare Visit Info ####
############################

# total number of combinations: \sum_i=1^k choose(n, k) + 1 for no cat
# so here it is sum(choose(4, 1:4)) +1

setting_labels <- expand.grid(outpatient = c("Outpatient", NA),
                    inpatient = c("Inpatient", NA),
                    ed = c("ED", NA),
                    obs_stay = c("Observational Stay", NA)) %>% 
  unite(., col = setting_label, sep = " + ", na.rm = T, remove =F) %>% 
  mutate(across(outpatient:obs_stay, ~ifelse(is.na(.), F, T))) %>% 
  mutate(setting_label = ifelse(setting_label == "", "Not any",
                                ifelse(setting_label == "Inpatient", "Inpatient Only",
                                       ifelse(setting_label == "Outpatient", "Outpatient Only",
                                              ifelse(setting_label == "ED", "ED Only",
                                                     ifelse(setting_label == "Observational Stay", "Observational Stay Only",
                                                            setting_label))))))
setting_labels <- setting_labels[order(rowSums(setting_labels[,2:5])), ]
setting_labels <- setting_labels %>% mutate(setting_label = factor(setting_label, 
                                 levels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")]),
                                 labels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")])))

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(patient_id,outpatient,ed,obs_stay,inpatient,admdate,oral_steroids, inhalers, antiacid_ppi:cough_suppressant)) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, oral_steroids, inhalers, antiacid_ppi:cough_suppressant) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(obs,patient_id,outpatient,ed,obs_stay,inpatient,admdate,oral_steroids, inhalers, antiacid_ppi:cough_suppressant)) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(obs,patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, oral_steroids, inhalers, antiacid_ppi:cough_suppressant) %>% 
  inner_join(setting_labels) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd <- sim_res_ssd %>% 
  inner_join(obs_locations)


full_reg_data <- nest(sim_res_ssd, .by = c("trial", "boot_trial"))
full_reg_data <- full_reg_data %>% 
  inner_join(boot_data %>% 
               mutate(boot_sample = map(boot_sample, ~inner_join(., index_locations))))

# all.equal(full_reg_data %>% slice(100) %>% select(boot_sample) %>% unnest(boot_sample),
#           full_reg_data %>% slice(1000) %>% select(boot_sample) %>% unnest(boot_sample))

print("dataset ready")

###########################
#### Regression Models ####
###########################

### Note exclude year 2001 and 2002 as no index cases occured then
# index_locations %>% count(year) %>% print(n = Inf)

#### Missed opportunities All --------------------------------------------------
# library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, setting_label, age_cat, female, rural, source, weekend, year, month,
           any_obesity)
  

  fit <- glm(miss ~ .,
             family = "binomial", data=reg_data)

  broom::tidy(fit)

  # fit <- logistf(miss ~ .,
  #                family = "binomial", data=reg_data,
  #                control = logistf.control( maxit = 1000),
  #                plcontrol = logistpl.control( maxit = 1000))
  # 
  # tibble(term = names(coef(fit)),
  #        estimate = coef(fit))
  
}

miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# rm(cluster)
gc()

print("mod 1 finished")

#### Missed opportunities All -inpatient only ind --------------------------------------------------

get_miss_res_inpatient_ind <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, inpatient, age_cat, female, rural, source, weekend, year, month,
           any_obesity)
  
  
  fit <- glm(miss ~ .,
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

miss_opp_res_inpatient_ind <- parallel::mclapply(1:max(full_reg_data$trial),
                                                  function(x){get_miss_res_inpatient_ind(x)}, 
                                                  mc.cores = num_cores)


miss_opp_res_inpatient_ind <- bind_rows(miss_opp_res_inpatient_ind) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# rm(cluster)
gc()

print("mod 2 finished")

#### Missed opportunities Medicaid ---------------------------------------------

get_miss_res_med <- function(trial_val){
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE))  %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    filter(source == "medicaid") %>% 
    select(miss, setting_label, age_cat, female, rural, race, weekend, year, month,
           any_obesity)
  
  
  fit <- glm(miss~.,
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
                                            function(x){get_miss_res_med(x)}, 
                                            mc.cores = num_cores)

miss_opp_res_med <- bind_rows(miss_opp_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

gc()

print("mod 3 finished")

#### Missed opportunities Medicaid - Inpatient_only ---------------------------------------------

get_miss_res_med_inpatient_ind <- function(trial_val){
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE))  %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    filter(source == "medicaid") %>% 
    select(miss, inpatient, age_cat, female, rural, race, weekend, year, month,
           any_obesity)
  
  
  fit <- glm(miss ~ .,
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res_med_inpatient_ind  <- parallel::mclapply(1:max(full_reg_data$trial),
                                       function(x){get_miss_res_med_inpatient_ind(x)}, 
                                       mc.cores = num_cores)

miss_opp_res_med_inpatient_ind <- bind_rows(miss_opp_res_med_inpatient_ind) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

gc()

print("mod 4 finished")

# compute duration by simulation (note this will be used in the next step)
full_reg_data_dur <- full_reg_data %>% 
  mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")) %>% 
                      group_by(patient_id) %>% 
                      summarise(duration = -min(days_since_index)) %>% 
                      ungroup()))


#### Delay duration All - LM --------------------------------------------------------

get_dur_res <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(duration, age_cat, female, rural, source, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- glm(duration~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                       function(x){get_dur_res(x)}, 
                                       mc.cores = num_cores)


miss_dur_res <- bind_rows(miss_dur_res) %>% 
  group_by(term) %>% 
  summarise(est = median(estimate, na.rm = T),
            low = quantile(estimate, probs = c(0.025), na.rm = T),
            high = quantile(estimate, probs = c(0.975), na.rm = T))


gc()

print("mod 5 finished")


#### Delay duration All - log-normal --------------------------------------------------------

get_dur_res_log_normal <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    mutate(log_duration = log(duration)) %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(log_duration, age_cat, female, rural, source, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- glm(log_duration~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_log_normal <- parallel::mclapply(1:max(full_reg_data$trial),
                                              function(x){get_dur_res_log_normal(x)}, 
                                              mc.cores = num_cores)


miss_dur_res_log_normal <- bind_rows(miss_dur_res_log_normal) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T))


gc()

print("mod 6 finished")


#### Delay duration All - weibull --------------------------------------------------------
library(survival)

get_dur_res_weibull <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    mutate(delta = 1) %>% 
    select(duration, delta, age_cat, female, rural, source, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- survreg(Surv(duration, delta) ~., dist = "weibull", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_weibull <- parallel::mclapply(1:max(full_reg_data$trial),
                                           function(x){get_dur_res_weibull(x)}, 
                                           mc.cores = num_cores)


miss_dur_res_weibull <- bind_rows(miss_dur_res_weibull) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T))


gc()

print("mod 7 finished")

#### Delay duration Medicaid ---------------------------------------------------

get_dur_res_med <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(source == "medicaid") %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(duration, age_cat, female, rural, race, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- glm(duration~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_dur_res_med(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_med <- bind_rows(miss_dur_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(estimate, na.rm = T),
            low = quantile(estimate ,probs = c(0.025), na.rm = T),
            high = quantile(estimate ,probs = c(0.975), na.rm = T))

gc()

print("mod 8 finished")


#### Delay duration Medicaid -log normal ---------------------------------------------------

get_dur_res_med_log_normal <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>%
    mutate(log_duration = log(duration)) %>% 
    filter(source == "medicaid") %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(log_duration, age_cat, female, rural, race, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- glm(log_duration~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_log_normal_med <- parallel::mclapply(1:max(full_reg_data$trial),
                                              function(x){get_dur_res_med_log_normal(x)}, 
                                              mc.cores = num_cores)


miss_dur_res_log_normal_med <- bind_rows(miss_dur_res_log_normal_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T))


gc()

print("mod 9 finished")


#### Delay duration Medicaid -weibull ---------------------------------------------------

get_dur_res_med_weibull <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>%
    filter(source == "medicaid") %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    mutate(delta = 1) %>% 
    select(duration, delta, age_cat, female, rural, race, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  fit <- survreg(Surv(duration, delta) ~., dist = "weibull", data=reg_data)
  
  broom::tidy(fit)
  
  
  
}


miss_dur_res_weibull_med <- parallel::mclapply(1:max(full_reg_data$trial),
                                                  function(x){get_dur_res_med_weibull(x)}, 
                                                  mc.cores = num_cores)


miss_dur_res_weibull_med <- bind_rows(miss_dur_res_weibull_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T))


gc()

print("mod 10 finished")

#### Delay Patient All ---------------------------------------------------------


get_delay_pat_res <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
    mutate(duration = replace_na(duration,0L)) %>%  
    mutate(miss = duration>0) %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  
  fit <- glm(miss~., family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_delay_pat_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                       function(x){get_delay_pat_res(x)}, 
                                       mc.cores = num_cores)

miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

gc()

print("mod 11 finished")

#### Delay Patient Medicaid ----------------------------------------------------


get_delay_pat_res_med <- function(trial_val){

  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
    mutate(duration = replace_na(duration,0L)) %>%  
    mutate(miss = duration>0) %>% 
    filter(source == "medicaid") %>% 
    filter(year>2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, race, oral_steroids_window, inhalers_window, year, month,
           any_obesity, antiacid_ppi_drugs_window:cough_suppressant_drugs_window)
  
  
  fit <- glm(miss~., family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

miss_delay_pat_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
                                         function(x){get_delay_pat_res_med(x)}, 
                                         mc.cores = num_cores)

miss_delay_pat_res_med <- bind_rows(miss_delay_pat_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


gc()

print("mod 12 finished")

### Alternative missed visits --------------------------------------------------
# 
# obs_locations2 <- tm %>% 
#   left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
#   filter(!is.na(obs)) %>% 
#   filter(setting_type !=4) %>%
#   distinct(patient_id,obs,setting_type,admdate) %>% 
#   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
#   mutate(dow=weekdays(as_date(admdate))) %>% 
#   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
#   select(patient_id,obs,setting,dow,weekend) 
# 
# index_locations2 <- tm %>% 
#   filter(days_since_index==0) %>% 
#   filter(setting_type !=4) %>%
#   distinct(patient_id,setting_type,admdate) %>% 
#   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
#   mutate(dow=weekdays(as_date(admdate))) %>% 
#   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
#   select(patient_id,setting,weekend,dow) 
# 
# sim_res_ssd <- sim_res_ssd %>% 
#   select(obs,trial) %>% 
#   inner_join(obs_locations2)
# 
# get_miss_res2 <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- sim_res_ssd %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% 
#                           mutate(miss=TRUE),
#                         index_locations2 %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo, by = "patient_id")
#   
#   
#   fit <- glm(miss~setting + age_cat + female + rural + source + weekend, family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# get_miss_res2(10)
# 
# cluster <- parallel::makeCluster(num_cores)
# # cluster <- parallel::makeCluster(6)
# 
# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("sim_res_ssd","get_miss_res2","index_locations2","reg_demo"),
#                         envir=environment())
# 
# 
# miss_opp_res2 <- parallel::parLapply(cl = cluster,
#                                     1:max(sim_res_ssd$trial),
#                                     function(x){get_miss_res2(x)})
# 
# 
# parallel::stopCluster(cluster)
# 
# miss_opp_res2 <- bind_rows(miss_opp_res2) %>% 
#   group_by(term) %>% 
#   summarise(est = mean(exp(estimate)),
#             low = quantile(exp(estimate),probs = c(0.025)),
#             high = quantile(exp(estimate),probs = c(0.975)))
# 
# 
# miss_opp_res2
# miss_opp_res

########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res,
                             miss_opp_res_inpatient_ind = miss_opp_res_inpatient_ind,
                             miss_opp_res_med = miss_opp_res_med,
                             miss_opp_res_med_inpatient_ind = miss_opp_res_med_inpatient_ind,
                             miss_dur_res = miss_dur_res,
                             miss_dur_res_log_normal = miss_dur_res_log_normal,
                             miss_dur_res_weibull = miss_dur_res_weibull,
                             miss_dur_res_med = miss_dur_res_med,
                             miss_dur_res_log_normal_med = miss_dur_res_log_normal_med,
                             miss_dur_res_weibull_med = miss_dur_res_weibull_med,
                             miss_delay_pat_res = miss_delay_pat_res,
                             miss_delay_pat_res_med = miss_delay_pat_res_med)

save(ssd_miss_risk_models, index_locations, file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/", proj_name, "/risk_model_report.Rmd"),
                  params = list(cond = tools::toTitleCase(proj_name)),
                  output_dir = out_path,
                  output_file = paste0(proj_name, "_risk_model_report_", lubridate::today() %>% format('%m-%d-%Y')))

