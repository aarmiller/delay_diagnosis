
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "pe"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/pe/",
                     base_path = "/Shared/Statepi_Diagnosis/prelim_results/pe/",  # base path to original prelim extract results
                     out_path = "/Shared/Statepi_Diagnosis/projects/pe/",   # path to output delay new results
                     ssd_name = "pe",
                     cp = 98,
                     upper_bound = 365,
                     final_model = "quad",
                     periodicity = TRUE,
                     boot_trials = 100,
                     sim_trials = 100,
                     miss_bins = c(1,2,3,4,5),
                     duration_bins = c(1,2,3,4,5,6,7,10,14,17,21))

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_base_path,"ssd_visit/")

out_path <- paste0(delay_base_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Connect to db
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

## Build the primary index_dx_visits -------------------------------------------
# Connect to dengue DB
db <- con

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()

# Set enrollment minimum
enroll_min <- delay_params$upper_bound

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

index_cases <- enrolid_crosswalk %>%
  inner_join(index_dx_dates) %>% 
  select(patient_id,enrolid,index_date) %>% 
  mutate(shift=0L)

num_cores <- 20
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

load(paste0(sim_out_path,"/sim_res.RData"))
sim_res_ssd <- sim_res
tmp1 <- sim_res_ssd %>% distinct(obs)

sim_obs_reduced <- sim_obs %>% 
  inner_join(tmp1, by = "obs") %>% 
  select(obs,patient_id,days_since_index)

sim_res_sim_obs <- sim_res_ssd %>% 
  inner_join(sim_obs_reduced, by = "obs") 

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


## Prepare antibiotic indicators -----------------------------------------------

rx_visits <- con %>% tbl("all_rx_visits") %>% 
  filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
  collect()

rx_visits <- rx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))

gc()

abx_codes <- read_csv(paste0("/Shared/Statepi_Diagnosis/projects/dengue/", "antibiotics.csv"))
abx_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% abx_codes$ndcnum) 

reg_demo <- reg_demo %>% 
  left_join(abx_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(abx_window=1L)) %>% 
  mutate_at(vars(abx_window),~replace_na(.,0L))

## Add abx to time map
abx_rx_vis <- abx_rx_vis %>% distinct(patient_id, admdate = date)

tm <- tm %>% rename(admdate = svcdate) %>% left_join(abx_rx_vis %>% mutate(abx = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(abx),~replace_na(.,0L)) 


## Prepare inhaler indicators --------------------------------------------------

# load(paste0(delay_params$out_path,"sim_results/reg_vars.RData"))
# reg_demo <- reg_demo %>%
#   select(-inhaler)

tmp <- read_csv("/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/data/all_inhaler_codes.csv")

inhaler_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% tmp$NDCNUM)

reg_demo <- reg_demo %>% 
  left_join(inhaler_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(inhaler_window=1L)) %>% 
  mutate_at(vars(inhaler_window),~replace_na(.,0L))

## Add inhaler to time map
inhaler_rx_vis <- inhaler_rx_vis %>% distinct(patient_id, admdate = date)

tm <- tm %>% left_join(inhaler_rx_vis %>% mutate(inhaler = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(inhaler),~replace_na(.,0L)) 

## Prepare prior asthma and COPD indicators ------------------------------------

#load asthma/copd codes
old_ssds <- read_rds(paste0("/Shared/Statepi_Diagnosis/params/ssd_codes/condition_specific/ssd_tb.RDS"))
old_ssds <- old_ssds[c("asthma", "copd")]
asthma_COPD_codes <- tibble(dx = old_ssds$asthma$icd9_codes, dx_ver= 9L, cond = "Asthma") %>% 
  bind_rows(tibble(dx = old_ssds$asthma$icd10_codes, dx_ver= 10L, cond = "Asthma")) %>% 
  bind_rows(tibble(dx = old_ssds$copd$icd9_codes, dx_ver= 9L, cond = "COPD")) %>% 
  bind_rows(tibble(dx = old_ssds$copd$icd10_codes, dx_ver= 10L, cond = "COPD")) %>% 
  distinct()

dx_visits <- con %>% tbl("all_dx_visits") %>% collect()

dx_visits <- dx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)

#Asthma
asthma <- dx_visits %>% 
  inner_join(asthma_COPD_codes %>% filter( cond == "Asthma") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(asthma_prior_cp = 1L) 

#COPD
copd <- dx_visits %>% 
  inner_join(asthma_COPD_codes %>% filter( cond == "COPD") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(copd_prior_cp = 1L) 

reg_demo <- reg_demo %>%
  left_join(asthma) %>% 
  left_join(copd) %>% 
  mutate_at(vars(asthma_prior_cp, copd_prior_cp),~replace_na(.,0L))

## Prepare prior chest imaging -------------------------------------------------
chest_ct <-c("71260", "71250", "71270")
chest_xray <- c("71010", "71015", "71020", "71021", "71022", "71023",
                "71030", "71034", "71035", "71101", "71111", "71120",
                "71045", "71046", "71047", "71048")

proc_visits <- con %>% tbl("all_proc_visits") %>% collect()

proc_visits <- proc_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)

chest_ct_visits <- proc_visits %>% 
  filter(proc %in% chest_ct)%>% 
  distinct(patient_id) %>% mutate(chest_ct_prior_cp = 1)

chest_xray_visits <- proc_visits %>% 
  filter(proc %in% chest_xray) %>% 
  distinct(patient_id) %>% mutate(chest_xray_prior_cp = 1)

reg_demo <- reg_demo %>%
  left_join(chest_ct_visits) %>% 
  left_join(chest_xray_visits) %>% 
  mutate_at(vars(chest_ct_prior_cp, chest_xray_prior_cp),~replace_na(.,0L))


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
  distinct(patient_id,outpatient,ed,obs_stay,inpatient,admdate,abx, inhaler) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, abx, inhaler) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct(obs,patient_id,outpatient,ed,obs_stay,inpatient,admdate,abx, inhaler) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(obs,patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, abx, inhaler) %>% 
  inner_join(setting_labels) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd <- sim_res_ssd %>% 
  inner_join(obs_locations)

###########################
#### Regression Models ####
###########################

#### Missed opportunities All - inpatient only --------------------------------------------------


get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id")
  
  
  fit <- glm(miss~ inpatient + age_cat + female + rural + source + weekend + abx + inhaler + asthma_prior_cp +
               copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                            function(x){get_miss_res(x)}, 
                                            mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# rm(cluster)
gc()

sim_res_ssd_temp <- sim_res_ssd %>% group_by(trial) %>% distinct(patient_id)
sim_res_ssd_temp <- sim_res_ssd_temp %>% ungroup() %>% 
  mutate(miss = 1L)
#### Delay Patient All ---------------------------------------------------------


get_delay_pat_res <- function(trial_val){
  
  tmp1 <- sim_res_ssd_temp %>% filter(trial==trial_val) 
  
  reg_data <- reg_demo %>% 
    left_join(tmp1, by = "patient_id") %>% 
    mutate(miss = replace_na(miss,0L))
  
  fit <- glm(miss~age_cat + female + rural + source + abx_window + inhaler_window + asthma_prior_cp +
               copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


# get_delay_pat_res(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("get_delay_pat_res"),
#                         envir=environment())
# 
# 
# miss_delay_pat_res <- parallel::parLapply(cl = cluster,
#                                           1:n_trials,
#                                           function(x){get_delay_pat_res(x)})

# parallel::stopCluster(cluster)

miss_delay_pat_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                       function(x){get_delay_pat_res(x)}, 
                                       mc.cores = num_cores)

miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

gc()

########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res,
                             miss_delay_pat_res = miss_delay_pat_res)



out_path <- paste0(delay_params$base_path,"risk_models_prelim/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

save(ssd_miss_risk_models,file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/pe/risk_model_report_prelim.Rmd",
                  params = list(cond = "Pulmonary Embolism"),
                  output_dir = out_path,
                  output_file = paste0(proj_name, "_risk_model_report_", lubridate::today() %>% format('%m-%d-%Y')))
