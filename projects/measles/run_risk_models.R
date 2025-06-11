
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "measles"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 14 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

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
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),"index_cases.RData"))

num_cores <- 50
cp_selected <- delay_params$cp - 1 # minus 1 as the risk factors are for within delay window

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

# age_cats <- c(-1,17,34,44,54,64,130)
age_cats <- c(-1,1,4,11,17,35,45,55,65,130) # pertussis age cats

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,age,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))


rx_visits <- con %>% tbl("all_rx_visits") %>% 
  filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
  collect()

rx_visits <- rx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))

gc()


tm <- tm %>% rename(admdate = svcdate)

## Prepare abx indicators --------------------------------------------------

# all_abx_codes <- read_csv(paste0(delay_params$out_path,"all_abx_codes.csv"))
all_abx_codes <- read_csv(paste0(stringr::str_replace(delay_params$out_path,
                                     paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
                "all_abx_codes.csv"))
all_abx_codes <- all_abx_codes %>%
  filter(!drug_name %in% tolower(c("Nitrofurantoin", "Metronidazole", "Trimethoprim",
                                   "Cefadroxil", "Rifampin", "Methenamine", "Demeclocycline",
                                   "Rifaximin", "Dapsone", "Gemifloxacin", 
                                   "Norfloxacin", "Fosfomycin", "Ofloxacin",
                                   "Carbenicillin", "Gatifloxacin", "Cefmetazole",
                                   "Cillin", "Fidaxomicin", "Sulfacetamide", "Sulfadiazine",
                                   "Sulfadimethoxine", "Sulfadoxine"))) %>% 
  distinct(ndcnum)

abx_rx_ids <- rx_visits %>% 
  filter(ndcnum %in% all_abx_codes$ndcnum) %>% 
  filter(date<index_date & date>=(index_date-cp_selected)) %>%
  distinct(patient_id) %>% 
  mutate(resp_antibiotic_drugs_window=1L)

reg_demo <- reg_demo %>%
  left_join(abx_rx_ids) %>% 
  mutate_at(vars(resp_antibiotic_drugs_window),~replace_na(.,0L))


## Prepare prior measles vaccination codes  ------------------------------------

# Hi Alan, 
# 
# Here are the different CPT codes for the measles vaccination. 
# 
# 90705 Measles virus vaccine, live, for subcutaneous use     
# 90706 Rubella virus vaccine, live, for subcutaneous use     
# 90707 Measles, mumps and rubella virus vaccine (MMR), live, for subcutaneous use   
# 90708 Measles and rubella virus vaccine, live, for subcutaneous use    
# 90710 Measles, mumps, rubella, and varicella vaccine (MMRV), live, for subcutaneous use     
# 
# Thanks,
# 
# Phil

vaccination <- c("90705", "90706", "90707", "90708", "90710")

proc_visits <- con %>% tbl("all_proc_visits") %>% collect()

proc_visits <- proc_visits %>%
  inner_join(select(demo1,patient_id,index_date)) %>%
  filter(days_since_index< (-cp_selected)) # any prior history of vaccination prior to cp

vaccination_hist <- proc_visits %>%
  filter(proc %in% vaccination) %>%
  distinct(patient_id) %>% mutate(vaccination_prior_cp = 1)

reg_demo <- reg_demo %>%
  left_join(vaccination_hist) %>%
  mutate_at(vars(vaccination_prior_cp),~replace_na(.,0L))


# #load asthma/copd codes
# old_ssds <- read_rds(paste0("/Shared/Statepi_Diagnosis/params/ssd_codes/condition_specific/ssd_tb.RDS"))
# old_ssds <- old_ssds[c("asthma", "copd")]
# asthma_COPD_codes <- tibble(dx = old_ssds$asthma$icd9_codes, dx_ver= 9L, cond = "Asthma") %>% 
#   bind_rows(tibble(dx = old_ssds$asthma$icd10_codes, dx_ver= 10L, cond = "Asthma")) %>% 
#   bind_rows(tibble(dx = old_ssds$copd$icd9_codes, dx_ver= 9L, cond = "COPD")) %>% 
#   bind_rows(tibble(dx = old_ssds$copd$icd10_codes, dx_ver= 10L, cond = "COPD")) %>% 
#   distinct()
# 
# dx_visits <- con %>% tbl("all_dx_visits") %>% collect()
# 
# dx_visits <- dx_visits %>% 
#   inner_join(select(demo1,patient_id,index_date)) %>% 
#   filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)
# 
# #Asthma
# asthma <- dx_visits %>% 
#   inner_join(asthma_COPD_codes %>% filter( cond == "Asthma") %>% select(-cond)) %>% 
#   distinct(patient_id) %>% mutate(asthma_prior_cp = 1L) 
# 
# #COPD
# copd <- dx_visits %>% 
#   inner_join(asthma_COPD_codes %>% filter( cond == "COPD") %>% select(-cond)) %>% 
#   distinct(patient_id) %>% mutate(copd_prior_cp = 1L) 
# 
# reg_demo <- reg_demo %>%
#   left_join(asthma) %>% 
#   left_join(copd) %>% 
#   mutate_at(vars(asthma_prior_cp, copd_prior_cp),~replace_na(.,0L))
# 
# ## Prepare prior chest imaging -------------------------------------------------
# chest_ct <-c("71260", "71250", "71270")
# chest_xray <- c("71010", "71015", "71020", "71021", "71022", "71023",
#                 "71030", "71034", "71035", "71101", "71111", "71120",
#                 "71045", "71046", "71047", "71048")
# 
# proc_visits <- con %>% tbl("all_proc_visits") %>% collect()
# 
# proc_visits <- proc_visits %>% 
#   inner_join(select(demo1,patient_id,index_date)) %>% 
#   filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)
# 
# chest_ct_visits <- proc_visits %>% 
#   filter(proc %in% chest_ct)%>% 
#   distinct(patient_id) %>% mutate(chest_ct_prior_cp = 1)
# 
# chest_xray_visits <- proc_visits %>% 
#   filter(proc %in% chest_xray) %>% 
#   distinct(patient_id) %>% mutate(chest_xray_prior_cp = 1)
# 
# reg_demo <- reg_demo %>%
#   left_join(chest_ct_visits) %>% 
#   left_join(chest_xray_visits) %>% 
#   mutate_at(vars(chest_ct_prior_cp, chest_xray_prior_cp),~replace_na(.,0L))
# 
# ## Prepare IPD 
# 
# #IPD
# ipd_codes <- codeBuildr::load_disease_codes("ipd")
# ipd_codes <- tibble(dx_ver = 9, dx = ipd_codes$ipd$icd9_codes) %>% 
#   bind_rows(tibble(dx_ver = 10, dx = ipd_codes$ipd$icd10_codes))
# 
# IPD <- dx_visits %>% inner_join(ipd_codes) %>%   
#   distinct(patient_id) %>% mutate(IPD_prior_cp = 1L) 
# 
# reg_demo <- reg_demo %>%
#   left_join(IPD) %>% 
#   mutate_at(vars(IPD_prior_cp),~replace_na(.,0L))

## Add in ID consults ----------------------------------------------------------
ID_consult_vis <- tm_stdprov %>% 
  filter(patient_id %in% local(index_cases$patient_id)) %>%
  filter(stdprov %in%  c("285", "448")) %>%
  distinct() %>%
  inner_join(index_cases %>%
               mutate(index_date = index_date + shift) %>%
               select(patient_id, index_date), by = "patient_id") %>%
  rename(admdate = svcdate) %>%
  filter(admdate<=index_date & (admdate)>=(index_date-(delay_params$upper_bound))) %>%
  distinct(patient_id, admdate)

tm <- tm  %>%
  left_join(ID_consult_vis %>% mutate(ID_consult = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(ID_consult),~replace_na(.,0L))


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
  distinct_at(vars(patient_id,outpatient,ed,obs_stay,inpatient,admdate, ID_consult)) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, ID_consult) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(obs,patient_id,outpatient,ed,obs_stay,inpatient,admdate, ID_consult)) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(obs,patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, ID_consult) %>% 
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

# exclude years 2001 and 2002 due to sparsity

#### Missed opportunities All (NO SETTING RISK FACTOR) --------------------------------------------------

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
    filter(as.numeric(year) > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, weekend, year, month, vaccination_prior_cp,
           ID_consult) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  fit <- glm(miss ~ .,
             family = "binomial", data=reg_data)
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  
}

miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

# rm(cluster)
gc()

print("mod 1a finished")

# #### Missed opportunities All --------------------------------------------------
# library(logistf)
# get_miss_res <- function(trial_val){
# 
#   # trial_val <- 1
# 
#   tmp1 <- full_reg_data %>%
#     filter(trial==trial_val)
# 
#   reg_data <- bind_rows(tmp1 %>% select(data) %>%
#                           unnest(data) %>%
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>%
#                           unnest(boot_sample) %>%
#                           mutate(miss=FALSE)) %>%
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(as.numeric(year) > 2002) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>%
#     select(miss, setting_label, age_cat, female, rural, source, weekend, year, month, vaccination_prior_cp,
#           ID_consult) %>%
#     drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
# 
#   fit <- logistf(miss ~ .,
#                  family = "binomial", data=reg_data,
#                  control = logistf.control( maxit = 1000),
#                  plcontrol = logistpl.control( maxit = 1000))
# 
#   tibble(term = names(coef(fit)),
#          estimate = coef(fit))
# 
# }
# 
# miss_opp_res_settings <- parallel::mclapply(1:max(full_reg_data$trial),
#                                    function(x){get_miss_res(x)},
#                                    mc.cores = num_cores)
# 
# 
# miss_opp_res_settings <- bind_rows(miss_opp_res_settings) %>%
#   group_by(term) %>%
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
#             low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
#             high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))
# 
# # rm(cluster)
# gc()
# 
# print("mod 1b finished")

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
    filter(as.numeric(year) > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>%
    select(miss, inpatient, age_cat, female, rural, source, weekend, year, month, vaccination_prior_cp,
           ID_consult) %>%
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data


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
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

# rm(cluster)
gc()

print("mod 1c finished")

#### Delay duration All --------------------------------------------------------

# compute duration by simulation (note this will be used in the next step)
full_reg_data_dur <- full_reg_data %>% 
  mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")) %>% 
                      group_by(patient_id) %>% 
                      summarise(duration = -min(days_since_index)) %>% 
                      ungroup()))

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
    filter(year > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(duration, age_cat, female, rural, source, year, month, vaccination_prior_cp,
           resp_antibiotic_drugs_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
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
            high = quantile(estimate, probs = c(0.975), na.rm = T),
            low_90 = quantile(estimate, probs = c(0.05), na.rm = T),
            high_90 = quantile(estimate, probs = c(0.95), na.rm = T) )


gc()

print("mod 2 finished")


# log normal model
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
    filter(year > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(log_duration, age_cat, female, rural, source, year, month, vaccination_prior_cp,
           resp_antibiotic_drugs_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
  fit <- glm(log_duration ~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_log_normal <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_dur_res_log_normal(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_log_normal <- bind_rows(miss_dur_res_log_normal) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))


gc()

print("mod 3 finished")

# weibul model
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
    filter(year > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    mutate(delta = 1) %>% 
    select(duration, delta, age_cat, female, rural, source, year, month, vaccination_prior_cp,
           resp_antibiotic_drugs_window) %>% 
    drop_na() %>%  # have to do at the end as obs not in bootsample and boot_id not in data
    mutate(source = as.factor(as.character(source))) # for some reason survreg will not drop the missing factor level medicaid
                                       # and when your run summary fit the NA coef will be given a value of 0
  
  fit <- survreg(Surv(duration, delta) ~., dist = "weibull", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_weibull <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_dur_res_weibull(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_weibull <- bind_rows(miss_dur_res_weibull) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

gc()

print("mod 4 finished")

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
    filter(year > 2002) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, year, month, vaccination_prior_cp,
           resp_antibiotic_drugs_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
  
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
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

gc()

print("mod 5 finished")


########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res,
                             # miss_opp_res_settings = miss_opp_res_settings,
                             miss_opp_res_inpatient_ind = miss_opp_res_inpatient_ind,
                             miss_dur_res = miss_dur_res,
                             miss_dur_res_log_normal = miss_dur_res_log_normal,
                             miss_dur_res_weibull = miss_dur_res_weibull,
                             miss_delay_pat_res = miss_delay_pat_res
                             )


save(ssd_miss_risk_models, index_locations, file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/", proj_name, "/risk_model_report.Rmd"),
                  params = list(cond = tools::toTitleCase(proj_name)),
                  output_dir = out_path,
                  output_file = paste0(proj_name, "_risk_model_report_", lubridate::today() %>% format('%m-%d-%Y')))

