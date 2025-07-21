
rm(list=ls()[!(ls() %in% c("proj_name", "project_name", "cps_index", "cps_vector", "delay_params_org"))])
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)
library(codeBuildr)

setwd("/Shared/Statepi_Diagnosis/atlan")

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

# args = commandArgs(trailingOnly=TRUE)

# name of condition
# proj_name <- args[1]
# proj_name <- "dengue"

load("/Shared/AML/params/final_delay_params_kaiser.RData")

delay_params <- final_delay_params[[proj_name]]

data_in_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

######################
#### Prepare Data ####
######################

#### Load Data  ----------------------------------------------------------------
# Load in simulation results
load(paste0(sim_out_path,"sim_res_ssd.RData"))
load(paste0(sim_out_path,"sim_res_all.RData"))
load(paste0(sim_out_path,"sim_obs_reduced.RData")) # comes from final_sim.R
load(paste0(delay_params$base_path,"/delay_results/all_dx_visits.RData"))
load(paste0(delay_params$base_path,"/delay_results/delay_tm.RData"))
load(paste0(delay_params$base_path,"/delay_results/sim_obs.RData"))# sim_obs
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
            "index_cases.RData"))
load(paste0(sim_out_path,"boot_data.RData"))

# update all_dx_visits
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

# Subset to project specific patient ids
index_dx_dates <- index_cases
patient_ids <- index_cases %>% distinct(patient_id)
n_patients <- nrow(patient_ids)

# update tm
tm <- tm %>%   
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>%
  mutate(days_since_index = days_since_index - shift) %>%
  select(-shift) %>%
  filter(days_since_index<=0) %>%
  select(patient_id, svcdate, days_since_index, outpatient, ed, inpatient, other) %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0)) # subset to only visits from AV, ED, IP, or IS

# update sim_obs
sim_obs <- sim_obs %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<0) %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index"))

sim_tm <- all_dx_visits %>%
  mutate(period = -days_since_index) %>%
  distinct(patient_id,period,days_since_index) %>%
  inner_join(sim_obs_reduced ,by = c("patient_id", "period"))

all_vis_count <- sim_tm %>%
  count(period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

# # merge observations into time map to extract visit types
obs_tm <- sim_tm %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,other,inpatient)

###################################
#### Analyze Bootstrap Results ####
###################################

compute_boot_stats <- function(sim_data,sim_obs_data,delay_params,n_patients){
  
  tmp1 <- sim_data %>% 
    inner_join(sim_obs_data, by = "obs") %>% 
    group_by(patient_id,boot_id) %>%
    summarise(duration=max(period),
              n_miss=n())
  
  count_duration <- function(val){
    tmp1 %>%
      ungroup() %>% 
      summarise(n=sum(duration>=val))
  }
  
  res1 <- tmp1  %>%
    ungroup() %>% 
    summarise(n_pat = n(),
              mean_dur = mean(duration),
              median_dur = median(duration),
              mean_n_miss = mean(n_miss),
              median_n_miss = median(n_miss))
  
  res2 <- tibble(duration_bin = delay_params$duration_bins) %>%
    mutate(res = map(duration_bin,count_duration)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp1),
           pct_all = 100*n/n_patients)
  
  count_miss <- function(val){
    tmp1 %>%
      ungroup() %>% 
      summarise(n = sum(n_miss>=val))
  }
  
  res3 <- tibble(miss_bin = delay_params$miss_bins) %>%
    mutate(res = map(miss_bin,count_miss)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp1),
           pct_all = 100*n/n_patients)
  
  return(list(stats = res1,
              dur_tab = res2,
              miss_tab = res3))
  
}

# setup cluster 
plan(multisession, workers = 10)

tmp1 <- sim_res_ssd %>% 
  # slice(1:30) %>% 
  mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
                                                         sim_obs_reduced,
                                                         delay_params = delay_params, 
                                                         n_patients = n_patients)))

tmp1 <- tmp1 %>% 
  select(sim_trial,boot_trial,sim_stats) %>% 
  mutate(stats = map(sim_stats,~.$stats)) %>% 
  mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>% 
  mutate(dur_tab = map(sim_stats,~.$dur_tab))

sim_stats_ssd <- tmp1 %>% 
  select(-sim_stats)

tmp2 <- sim_res_all %>% 
  # slice(1:30) %>% 
  mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
                                                         sim_obs_reduced,
                                                         delay_params = delay_params, 
                                                         n_patients = n_patients)))

tmp2 <- tmp2 %>% 
  select(sim_trial,boot_trial,sim_stats) %>% 
  mutate(stats = map(sim_stats,~.$stats)) %>% 
  mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>% 
  mutate(dur_tab = map(sim_stats,~.$dur_tab))

sim_stats_all <- tmp2 %>% 
  select(-sim_stats)

save(sim_stats_ssd,sim_stats_all,file = paste0(sim_out_path,"sim_stats.RData"))

plan(sequential)

###################################
#### Export Simulation Results ####
###################################

# Note: We may want to move parts of this section to a separate markdown report script

### Generate aggregate statistics ----------------------------------------------

aggregate_sim_stats <- function(sim_stats_data){
  tmp1 <- sim_stats_data %>% 
    select(sim_trial,boot_trial,stats) %>% 
    unnest(stats) %>% 
    summarise_at(vars(n_pat:median_n_miss),~round(mean(.),2))
  
  tmp2 <- sim_stats_data %>% 
    select(sim_trial,boot_trial,stats) %>% 
    unnest(stats) %>% 
    summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.025)),2))
  
  tmp3 <- sim_stats_data %>% 
    select(sim_trial,boot_trial,stats) %>% 
    unnest(stats) %>% 
    summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.975)),2))
  
  
  res1 <- inner_join(gather(tmp1, key=measure, value = mean),
                     gather(tmp2, key=measure, value = low)) %>% 
    inner_join(gather(tmp3, key=measure, value = high)) %>% 
    mutate(measure_out=paste0(mean," (",low,"-",high,")"))
  
  return(list(main_stats = res1))
}

# generate aggregated results of:
#   n_pat - number of patient missed
#   mean_dur - mean duration of misses
#   median_dur - median duration of misses
#   mean_n_miss - mean number of misses per patient
#   median_n_miss - median number of misses per patient
agg_stats_ssd <- aggregate_sim_stats(sim_stats_ssd) 
agg_stats_all <- aggregate_sim_stats(sim_stats_all)

### Generate bins for number of misses per patient -----------------------------

# SSD Visits
tmp1 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
  gather(key = key, value = mean, -miss_bin) %>% 
  mutate(mean = round(mean,2))

tmp2 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
  gather(key = key, value = low, -miss_bin) %>% 
  mutate(low = round(low,2))

tmp3 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
  gather(key = key, value = high, -miss_bin) %>% 
  mutate(high = round(high,2))

miss_bins_ssd <- inner_join(tmp1,tmp2) %>% 
  inner_join(tmp3) %>% 
  mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
  select(miss_bin,key,out) %>% 
  spread(key = key, value = out)

# All visits
tmp1 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
  gather(key = key, value = mean, -miss_bin) %>% 
  mutate(mean = round(mean,2))

tmp2 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
  gather(key = key, value = low, -miss_bin) %>% 
  mutate(low = round(low,2))

tmp3 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
  gather(key = key, value = high, -miss_bin) %>% 
  mutate(high = round(high,2))

miss_bins_all <- inner_join(tmp1,tmp2) %>% 
  inner_join(tmp3) %>% 
  mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
  select(miss_bin,key,out) %>% 
  spread(key = key, value = out)


### Generate bins for duration of misses ---------------------------------------

# SSD Visits
tmp1 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
  gather(key = key, value = mean, -duration_bin) %>% 
  mutate(mean = round(mean,2))

tmp2 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
  gather(key = key, value = low, -duration_bin) %>% 
  mutate(low = round(low,2))

tmp3 <- sim_stats_ssd %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
  gather(key = key, value = high, -duration_bin) %>% 
  mutate(high = round(high,2))

dur_bins_ssd <- inner_join(tmp1,tmp2) %>% 
  inner_join(tmp3) %>% 
  mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
  select(duration_bin,key,out) %>% 
  spread(key = key, value = out)

# ALL Visits
tmp1 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
  gather(key = key, value = mean, -duration_bin) %>% 
  mutate(mean = round(mean,2))

tmp2 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
  gather(key = key, value = low, -duration_bin) %>% 
  mutate(low = round(low,2))

tmp3 <- sim_stats_all %>% 
  select(sim_trial,boot_trial,dur_tab) %>% 
  unnest(dur_tab) %>% 
  group_by(duration_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
  gather(key = key, value = high, -duration_bin) %>% 
  mutate(high = round(high,2))

dur_bins_all <- inner_join(tmp1,tmp2) %>% 
  inner_join(tmp3) %>% 
  mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
  select(duration_bin,key,out) %>% 
  spread(key = key, value = out)

### Setting summary miss opportunities --------------------------------------------

# Load in the setting types
load(paste0(data_in_path,"visit_info.RData"))

# update tm_stdplac
tm_stdplac <- tm_stdplac %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

tm_stdplac <- tm_stdplac %>% 
  left_join(select(sim_tm,-period))

encounter_subtype_labels <- tribble(~encounter_subtype,~label, 
                                    "AI", "Acute inpatient stay",
                                    "CH", "Online Text Chat",
                                    "DI"," Dialysis",
                                    "EM", "Email/Patient Message",
                                    "HA", "Hospital ambulatory",
                                    "HH", "Home health",
                                    "HS", "Hospice",
                                    "NH", "Nursing home",
                                    "OA", "Other Asynchronous Virtual Care",
                                    "OI", "Online Intake/Triage Form",
                                    "OB", "Observation bed",
                                    "OC", "Outpatient clinic visit",
                                    "PP", "Provider to Provider/Facility to Facility Virtual Care",
                                    "RH", "Rehab",
                                    "SD", "Same-day surgery",
                                    "SN", "Skilled nursing facility",
                                    "TC", "Telephone - Call Center/Advice/Triage",
                                    "TN", "Telephone - Not Specified",
                                    "TS", "Telephone - Scheduled Visit",
                                    "VV", "Video Visit",
                                    "UC", "Urgent care",
                                    "OT", "Other non-hospital")


generate_setting_counts_2_kaiser <- function(obs_tm, bootstrap_data, sim_res, tm_stdplac){
  
  tmp <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>% 
    inner_join(distinct(tm_stdplac,obs,encounter_subtype),by = join_by(obs),relationship = "many-to-many")
  
  tmp2 <- tmp %>%
    group_by(trial) %>%
    count(encounter_subtype) %>%
    ungroup()
  
  tmp_fill <- tmp2 %>%
    distinct(encounter_subtype) %>%
    mutate(trial = map(encounter_subtype,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
    unnest(trial)
  
  stdplac_res <- tmp2 %>%
    left_join(tmp_fill,.,by = c("encounter_subtype", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    ungroup(trial) %>% 
    inner_join(sim_res$sim_res_obs %>% distinct(boot_trial, trial) %>% 
                 inner_join(bootstrap_data, by = "boot_trial") %>% 
                 distinct(trial, n_miss), by = "trial") %>% 
    mutate(pct_miss_vis_days = n/n_miss*100) %>% 
    group_by(encounter_subtype) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975),
              pct_mean_vis_days = mean(pct_miss_vis_days),
              pct_vis_days_low = quantile(pct_miss_vis_days,probs = 0.025),
              pct_vis_days_high = quantile(pct_miss_vis_days,probs = 0.975)) %>%
    # mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
    #        pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    # select(stdplac,n,pct,everything())
    select(encounter_subtype, everything()) %>% 
    left_join(encounter_subtype_labels, by = "encounter_subtype") 
  
  
  tmp2 <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>%
    group_by(trial) %>%
    summarise(outpatient = sum(outpatient),
              ed = sum(ed),
              other = sum(other),
              inpatient = sum(inpatient)) %>% 
    ungroup() %>% 
    gather(key = setting_type, value=n, -trial)
  
  tmp_fill <- tmp2 %>%
    distinct(setting_type) %>%
    mutate(trial = map(setting_type,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
    unnest(trial)
  
  setting_type_res <- tmp2 %>%
    left_join(tmp_fill,.,by = c("setting_type", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    ungroup(trial) %>% 
    inner_join(sim_res$sim_res_obs %>% distinct(boot_trial, trial) %>% 
                 inner_join(bootstrap_data, by = "boot_trial") %>% 
                 distinct(trial, n_miss), by = "trial") %>% 
    mutate(pct_miss_vis_days = n/n_miss*100) %>% 
    group_by(setting_type) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975),
              pct_mean_vis_days = mean(pct_miss_vis_days),
              pct_vis_days_low = quantile(pct_miss_vis_days,probs = 0.025),
              pct_vis_days_high = quantile(pct_miss_vis_days,probs = 0.975)) %>%
    # mutate(setting_type = setting_type_labels(setting_type)) %>%
    # mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
    #        pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    # select(setting_type,n,pct,everything())
    select(setting_type, everything())
  
  setting_type_res_p2 <- tmp2 %>%
    left_join(tmp_fill,.,by = c("setting_type", "trial")) %>%
    mutate(n=replace_na(n,0)) %>% 
    inner_join(sim_res$sim_res_obs %>% distinct(boot_trial, trial) %>% 
                 inner_join(bootstrap_data, by = "boot_trial") %>% 
                 distinct(trial, potential_miss_setting) %>% 
                 unnest(potential_miss_setting), by = c("trial", "setting_type")) %>% 
    mutate(pct_potential_miss_vis_days_setting = n/n_potential_miss*100) %>% 
    group_by(setting_type) %>%
    summarise(n_mean_potential_miss_vis_days_setting = mean(n_potential_miss),
              n_potential_miss_vis_days_setting_low = quantile(n_potential_miss,probs = 0.025),
              n_potential_miss_vis_days_setting_high = quantile(n_potential_miss,probs = 0.975),
              pct_mean_potential_miss_vis_days_setting = mean(pct_potential_miss_vis_days_setting),
              pct_potential_miss_vis_days_setting_low = quantile(pct_potential_miss_vis_days_setting,probs = 0.025),
              pct_potential_miss_vis_days_setting_high = quantile(pct_potential_miss_vis_days_setting,probs = 0.975)) 
  
  setting_type_res <- setting_type_res %>% inner_join(setting_type_res_p2)
  
  return(list(stdplac_res = stdplac_res,
              setting_type_res = setting_type_res))
}


sim_res <- list(sim_res_obs= sim_res_ssd %>% mutate(trial = row_number()) %>% 
                  unnest(cols = c(res)),
                sim_res_miss_tab = sim_res_ssd %>% mutate(trial = row_number()) %>%
                  distinct(trial))

setting_counts_ssd <- generate_setting_counts_2_kaiser(obs_tm = obs_tm,
                                                bootstrap_data = boot_data %>% select(boot_trial, n_miss = n_miss_ssd, potential_miss_setting = potential_miss_setting_ssd),
                                                sim_res = sim_res,
                                                tm_stdplac = tm_stdplac)

sim_res <- list(sim_res_obs= sim_res_all %>% mutate(trial = row_number()) %>% 
                  unnest(cols = c(res)),
                sim_res_miss_tab = sim_res_all %>% mutate(trial = row_number()) %>%
                  distinct(trial))

setting_counts_all <- generate_setting_counts_2_kaiser(obs_tm = obs_tm,
                                                bootstrap_data = boot_data %>% select(boot_trial, n_miss = n_miss_all, potential_miss_setting = potential_miss_setting_all),
                                                sim_res = sim_res,
                                                tm_stdplac = tm_stdplac)

### Setting summary index ------------------------------------------------------


# generate counts by stdplac version 2
generate_setting_counts_2_kaiser <- function(tm_data,sim_res_data, bootstrap_data, sim_res_sim_obs_data, tm_stdplac_data){
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()
  
  index_stdplac_counts_temp <- tm_data %>% 
    filter(days_since_index==0) %>%
    inner_join(tm_stdplac_data,by = join_by(patient_id, svcdate, days_since_index)) %>% 
    # filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(encounter_subtype_labels, by = "encounter_subtype") %>% 
    distinct(patient_id,svcdate,label) %>% 
    inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
    group_by(boot_trial) %>% 
    count(label,name = "index_visits") %>% 
    ungroup()
  
  tmp_fill <- index_stdplac_counts_temp %>%
    distinct(label) %>%
    mutate(boot_trial = map(label,~1:max(index_stdplac_counts_temp$boot_trial))) %>%
    unnest(boot_trial)
  
  index_stdplac_counts <- index_stdplac_counts_temp %>% 
    left_join(tmp_fill,.,by = c("label", "boot_trial")) %>%
    mutate(index_visits=replace_na(index_visits,0)) %>%
    group_by(boot_trial) %>% 
    mutate(index_pct1_temp = 100*index_visits/sum(index_visits),
           index_pct2_temp = 100*index_visits/total_patients) %>% 
    ungroup() %>% 
    group_by(label) %>% 
    summarise(index_count_mean = mean(index_visits),
              index_count_lowerCI = quantile(index_visits, probs = 0.025),
              index_count_upperCI = quantile(index_visits, probs = 0.975),
              index_pct1_mean = mean(index_pct1_temp),   # percent of total index locations,
              index_pct1_lowerCI = quantile(index_pct1_temp, probs = 0.025), # percent of total index locations lower CI
              index_pct1_upperCI = quantile(index_pct1_temp, probs = 0.975), # percent of total index locations upper CI
              index_pct2_mean = mean(index_pct2_temp), # percent of total patients %>% # percent of total index locations w/o Other
              index_pct2_lowerCI = quantile(index_pct2_temp, probs = 0.025), # percent of total patients %>% # percent of total index locations w/o Other lower CI
              index_pct2_upperCI = quantile(index_pct2_temp, probs = 0.975)) # percent of total patients %>% # percent of total index locations w/o Other lower CI%>% 
  
  ## get observation timemap -----------------------------------------------------
  tmp  <- sim_res_data %>% 
    inner_join(sim_res_sim_obs_data, by = "obs") %>% 
    inner_join(tm_stdplac_data %>% filter(!is.na(obs)) %>% select(obs, encounter_subtype), by = "obs",relationship = "many-to-many") %>% # na obs reflects non select sim obs
    # filter(!(stdplac %in% c(81,41,42))) %>% 
    left_join(encounter_subtype_labels, by = "encounter_subtype") 
  
  tmp_fill <- tmp %>%
    distinct(label) %>%
    mutate(trial = map(label,~1:max(sim_res_data$trial))) %>%
    unnest(trial)
  
  sim_res_setting_counts <- tmp %>% 
    group_by(trial, label) %>% 
    summarise(n=n(),
              duration = mean(-days_since_index)) %>%  
    ungroup() %>% 
    left_join(tmp_fill,.,by = c("label", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    full_join(sim_res_data %>% distinct(trial, boot_trial) %>% 
                inner_join(index_stdplac_counts_temp %>% mutate(n_index = index_visits), by = "boot_trial") %>% 
                distinct(trial, label, n_index), by = c("trial", "label")) %>% 
    mutate(n_index=replace_na(n_index,0)) %>%
    mutate(n=replace_na(n,0)) %>%
    mutate(pct_opp=replace_na(pct_opp,0)) %>%
    mutate(total_opps = n+n_index) %>% 
    mutate(pct_opp_missed = 100*n/(total_opps)) %>% 
    group_by(label) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              dur_mean = mean(duration, na.rm = T),
              dur_low = quantile(duration, na.rm = T, probs = 0.025),
              dur_high = quantile(duration, na.rm = T, probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              total_opps_mean = mean(total_opps),
              total_opps_low = quantile(total_opps,probs = 0.025),
              total_opps_high = quantile(total_opps,probs = 0.975),
              pct_opp_missed_mean = mean(pct_opp_missed, na.rm = T),
              pct_opp_missed_low = quantile(pct_opp_missed, na.rm = T, probs = 0.025),
              pct_opp_missed_high = quantile(pct_opp_missed, na.rm = T, probs = 0.975)) 
  
  sim_res_setting_counts %>% 
    full_join(index_stdplac_counts,by = "label")
  
}


# generate counts by setting
generate_setting_counts_3_kaiser <- function(tm_data, sim_res_data, bootstrap_data, sim_res_sim_obs_data, tm_stdplac_data){
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()
  
  index_setting_counts_temp <- tm_data %>% 
    filter(days_since_index==0) %>%
    inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
    group_by(boot_trial) %>% 
    summarise(outpatient = sum(outpatient),
              ed = sum(ed),
              other= sum(other),
              inpatient = sum(inpatient)) %>% 
    ungroup() %>% 
    gather(key = Setting, value = n, -boot_trial)
  
  tmp_fill <- index_setting_counts_temp %>%
    distinct(Setting) %>%
    mutate(boot_trial = map(Setting,~1:max(index_setting_counts_temp$boot_trial))) %>%
    unnest(boot_trial)
  
  index_setting_counts <- index_setting_counts_temp %>% 
    left_join(tmp_fill,.,by = c("Setting", "boot_trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(boot_trial) %>% 
    mutate(index_pct1_temp = 100*n/sum(n),
           index_pct2_temp = 100*n/total_patients) %>% 
    ungroup() %>% 
    group_by(Setting) %>% 
    summarise(index_count_mean = mean(n),
              index_count_lowerCI = quantile(n, probs = 0.025),
              index_count_upperCI = quantile(n, probs = 0.975),
              index_pct1_mean = mean(index_pct1_temp),   # percent of total index locations,
              index_pct1_lowerCI = quantile(index_pct1_temp, probs = 0.025), # percent of total index locations lower CI
              index_pct1_upperCI = quantile(index_pct1_temp, probs = 0.975), # percent of total index locations upper CI
              index_pct2_mean = mean(index_pct2_temp), # percent of total patients %>% # percent of total index locations w/o Other
              index_pct2_lowerCI = quantile(index_pct2_temp, probs = 0.025), # percent of total patients %>% # percent of total index locations w/o Other lower CI
              index_pct2_upperCI = quantile(index_pct2_temp, probs = 0.975)) # percent of total patients %>% # percent of total index locations w/o Other lower CI%>% 
  
  ## get observation timemap -----------------------------------------------------
  tmp2  <- sim_res_data %>% 
    inner_join(sim_res_sim_obs_data, by = "obs") %>% 
    inner_join(tm_data, by = join_by(patient_id, days_since_index)) %>% 
    select(trial, days_since_index, outpatient, ed, other, inpatient)
  
  tmp_fill2 <- tibble(Setting = names(tmp2)[-(1:2)]) %>% 
    mutate(trial = map(Setting,~1:max(sim_res_data$trial))) %>%
    unnest(trial)
  
  sim_res_setting_counts2 <- tibble(Setting = names(tmp2)[-(1:2)]) %>% 
    mutate(out = map(Setting, ~tmp2 %>% filter(get(.x) == 1) %>% 
                       group_by(trial) %>% 
                       summarise(n = sum(get(.x)),
                                 duration = mean(-days_since_index)) %>% 
                       ungroup())) %>% 
    unnest(out) %>% 
    left_join(tmp_fill2,.,by = c("Setting", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    full_join(sim_res_data %>% distinct(trial, boot_trial) %>% 
                inner_join(index_setting_counts_temp %>% mutate(n_index = n), by = "boot_trial") %>% 
                distinct(trial, Setting, n_index), by = c("trial", "Setting")) %>% 
    mutate(n_index=replace_na(n_index,0)) %>%
    mutate(n=replace_na(n,0)) %>%
    mutate(pct_opp=replace_na(pct_opp,0)) %>%
    mutate(total_opps = n+n_index) %>% 
    mutate(pct_opp_missed = 100*n/(total_opps)) %>% 
    group_by(Setting) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              dur_mean = mean(duration, na.rm = T),
              dur_low = quantile(duration, na.rm = T, probs = 0.025),
              dur_high = quantile(duration, na.rm = T, probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              total_opps_mean = mean(total_opps),
              total_opps_low = quantile(total_opps,probs = 0.025),
              total_opps_high = quantile(total_opps,probs = 0.975),
              pct_opp_missed_mean = mean(pct_opp_missed, na.rm = T),
              pct_opp_missed_low = quantile(pct_opp_missed, na.rm = T, probs = 0.025),
              pct_opp_missed_high = quantile(pct_opp_missed, na.rm = T, probs = 0.975)) 
  
  sim_res_setting_counts2 %>% 
    full_join(index_setting_counts,by = "Setting")
  
}

sim_res <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  unnest(cols = c(res))

# Generate the index counts by stdplac
setting_counts_index_by_stdplac <- generate_setting_counts_2_kaiser(tm_data = tm,
                                                  bootstrap_data = boot_data %>% select(boot_trial, boot_sample),
                                                  sim_res_data = sim_res,
                                                  sim_res_sim_obs_data = sim_obs_reduced %>% mutate(days_since_index = -period),
                                                  tm_stdplac_data = tm_stdplac)

# Generate the index counts by setting
setting_counts_index_by_setting <- generate_setting_counts_3_kaiser(tm_data = tm,
                                                             bootstrap_data = boot_data %>% select(boot_trial, boot_sample),
                                                             sim_res_data = sim_res,
                                                             sim_res_sim_obs_data = sim_obs_reduced %>% mutate(days_since_index = -period),
                                                             tm_stdplac_data = tm_stdplac)

### New setting counts_ssd table -----------------------------------------------

# define change_point
cp <- delay_params$cp

#load ssd codes
ssd_codes <- codeBuildr::load_ssd_codes(delay_params$ssd_name) %>% 
  filter(type %in% c("icd9","icd10")) %>% 
  mutate(dx_ver = ifelse(type=="icd9","09","10")) %>% 
  select(dx = code,dx_ver)

### Compute Index Setting Counts

out_index_n <- tm %>% 
  filter(days_since_index==0) %>% 
  distinct(patient_id, outpatient,ed,inpatient,other) %>% 
  gather(key = Setting, value = value, -patient_id) %>% 
  group_by(Setting) %>% 
  summarise(index_n = sum(value)) %>% 
  mutate(percent_index = index_n/sum(index_n)*100)

# out_index_n <- bind_rows(out_index_n,
#                          filter(out_index_n,Setting=="inpatient") %>% 
#                            mutate(Setting = "inpatient visit"))

### Compute Potential Opportunity Counts

# Count for visit days
out_pot_opps <- all_dx_visits %>% 
  inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  filter(between(days_since_index,-cp+1,-1)) %>% 
  distinct(obs,outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value, -obs) %>%
  group_by(Setting) %>% 
  summarise(potential_opps = sum(value))

# build counts for inpatient stays NO CASE ID

# tmp <- all_dx_visits %>% 
#   inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
#   distinct(patient_id,days_since_index) %>% 
#   inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
#   # inner_join(distinct(caseids,obs,caseid,patient_id),by = join_by(patient_id, obs)) %>% no caseid
#   filter(between(days_since_index,-cp+1,-1)) %>% 
#   # distinct(patient_id,caseid) %>%  no caseid
#   distinct(patient_id) %>% 
#   count(name = "potential_opps") %>% 
#   mutate(Setting = "inpatient visit")

# combine daily settings and inpatient visits
# out_pot_opps <- bind_rows(out_pot_opps,tmp)


### Compute Setting Miss Counts 

obs_tm <- sim_obs %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,other,inpatient)

# # merge observations into time map to extract visit types
tmp <- sim_res_ssd %>% 
  mutate(res = map(res,~inner_join(.,obs_tm,by = "obs")))

# get setting counts for each trial
setting_counts1 <- tmp %>% 
  mutate(n = map_int(res,nrow)) %>% 
  mutate(outpatient = map_int(res,~sum(.$outpatient))) %>% 
  mutate(ed = map_int(res,~sum(.$ed))) %>% 
  mutate(other = map_int(res,~sum(.$other))) %>% 
  mutate(inpatient = map_int(res,~sum(.$inpatient))) %>% 
  select(sim_trial,boot_trial,n:inpatient)

# compute miss counts for outpatient, ed, inpatient days and other
out1 <- setting_counts1 %>% 
  select(-n) %>% 
  gather(key = Setting, value = n, -sim_trial,-boot_trial) %>% 
  group_by(Setting) %>% 
  summarise(miss_mean = ceiling(mean(n)),
            miss_median = ceiling(median(n)),
            miss_low = ceiling(quantile(n,probs = 0.025)),
            miss_high = ceiling(quantile(n,probs = 0.975)))

# compute miss counts for inpatient visits NO CASEID

# # merge caseids into simulation results
# tmp <- sim_res_ssd %>% 
#   unnest(res) %>% 
#   inner_join(distinct(caseids,obs,caseid,patient_id), by = "obs") %>% 
#   distinct(sim_trial,boot_trial,patient_id,boot_id,caseid)
# 
# # compute counts for inpatient stays
# out2 <- tmp %>% 
#   count(sim_trial,boot_trial) %>% 
#   summarise(miss_mean = ceiling(mean(n)),
#             miss_median = ceiling(median(n)),
#             miss_low = ceiling(quantile(n,probs = 0.025)),
#             miss_high = ceiling(quantile(n,probs = 0.975))) %>% 
#   mutate(`Setting`= "inpatient visit")
# 
# out_miss_n <- bind_rows(out1,out2)

out_miss_n <- out1

# Compute Percentage Missed Visit Days (i.e. of all % of all missed visit days 
# during the delay window [i.e., observed - expected] days deemed to be missed 
# opportunities [here the denominator is the same across settings and since patients 
# can encounter more than one setting in a given day, these %’s will sum up to >= 100%]): 

out_miss_frac <- setting_counts1 %>% 
  mutate_at(vars(outpatient,ed,other,inpatient),~100*./n) %>% 
  select(outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean = mean(value),
            frac_median = median(value),
            frac_low = quantile(value,probs = 0.025),
            frac_high = quantile(value,probs = 0.975))


# Compute % of all missed opportunities by setting (here the %’s will sum to 100) 
out_miss_frac_miss_opps <- setting_counts1 %>% 
  select(-n) %>% 
  mutate(n = outpatient+ed+other+inpatient) %>% 
  mutate_at(vars(outpatient,ed,other,inpatient),~100*./n) %>% 
  select(outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean_all_miss_opp = mean(value),
            frac_median_all_miss_opp = median(value),
            frac_low_all_miss_opp = quantile(value,probs = 0.025),
            frac_high_all_miss_opp = quantile(value,probs = 0.975))


### Assemble table
# Final Setting Count Table
setting_table_new_raw_ssd <- tibble(Setting = c("outpatient","ed","other","inpatient")) %>% 
  left_join(out_index_n, by = "Setting") %>% 
  left_join(out_pot_opps, by = "Setting") %>% 
  left_join(out_miss_n, by = "Setting") %>% 
  left_join(out_miss_frac, by = "Setting") %>% 
  left_join(out_miss_frac_miss_opps, by = "Setting")

# Output Table for Reporting
setting_counts_index_by_setting_ssd_new <- setting_table_new_raw_ssd %>% 
  mutate(`% of all Index Locations` = round(percent_index, 2)) %>% 
  mutate(`Number of Missed Opportunities` = paste0(miss_mean, " (", miss_low, "-", miss_high, ")")) %>% 
  mutate(`% of Missed opportunities` = paste0(round(frac_mean_all_miss_opp,2), " (", round(frac_low_all_miss_opp,2), "-", round(frac_high_all_miss_opp,2), ")")) %>% 
  mutate(`% of Missed Visit Days` = paste0(round(frac_mean,2), " (", round(frac_low,2), "-", round(frac_high,2), ")")) %>% 
  select(Setting,
         `Index Count`=index_n, 
         `% of all Index Locations`,
         `Potential Missed Opportunity Visit Days`=potential_opps,
         `Number of Missed Opportunities`, 
         `% of Missed opportunities`,
         `% of Missed Visit Days`) %>% 
  mutate(`% of Missed Visit Days` = ifelse(Setting == "inpatient visit", "", `% of Missed Visit Days`),
         `% of Missed opportunities` = ifelse(Setting == "inpatient visit", "", `% of Missed opportunities`))


### New setting counts_all_visits table -----------------------------------------------

### Compute Potential Opportunity Counts

# Count for visit days
out_pot_opps <- all_dx_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  filter(between(days_since_index,-cp+1,-1)) %>% 
  distinct(obs,outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value, -obs) %>%
  group_by(Setting) %>% 
  summarise(potential_opps = sum(value))

# # build counts for inpatient stays NO CASE ID
# tmp <- all_dx_visits %>% 
#   distinct(patient_id,days_since_index) %>% 
#   inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
#   inner_join(distinct(caseids,obs,caseid,patient_id),by = join_by(patient_id, obs)) %>% 
#   filter(between(days_since_index,-cp+1,-1)) %>% 
#   distinct(patient_id,caseid) %>% 
#   count(name = "potential_opps") %>% 
#   mutate(Setting = "inpatient visit")
# 
# # combine daily settings and inpatient visits
# out_pot_opps <- bind_rows(out_pot_opps,tmp)


### Compute Setting Miss Counts 

obs_tm <- sim_obs %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,other,inpatient)

# # merge observations into time map to extract visit types
tmp <- sim_res_all %>% 
  mutate(res = map(res,~inner_join(.,obs_tm,by = "obs")))

# get setting counts for each trial
setting_counts1 <- tmp %>% 
  mutate(n = map_int(res,nrow)) %>% 
  mutate(outpatient = map_int(res,~sum(.$outpatient))) %>% 
  mutate(ed = map_int(res,~sum(.$ed))) %>% 
  mutate(other = map_int(res,~sum(.$other))) %>% 
  mutate(inpatient = map_int(res,~sum(.$inpatient))) %>% 
  select(sim_trial,boot_trial,n:inpatient)

# compute miss counts for outpatient, ed, inpatient days and other
out1 <- setting_counts1 %>% 
  select(-n) %>% 
  gather(key = Setting, value = n, -sim_trial,-boot_trial) %>% 
  group_by(Setting) %>% 
  summarise(miss_mean = ceiling(mean(n)),
            miss_median = ceiling(median(n)),
            miss_low = ceiling(quantile(n,probs = 0.025)),
            miss_high = ceiling(quantile(n,probs = 0.975)))

# compute miss counts for inpatient visits
# NO CASEID
# # merge caseids into simulation results
# tmp <- sim_res_all %>% 
#   unnest(res) %>% 
#   inner_join(distinct(caseids,obs,caseid,patient_id), by = "obs") %>% 
#   distinct(sim_trial,boot_trial,patient_id,boot_id,caseid)
# 
# # compute counts for inpatient stays
# out2 <- tmp %>% 
#   count(sim_trial,boot_trial) %>% 
#   summarise(miss_mean = ceiling(mean(n)),
#             miss_median = ceiling(median(n)),
#             miss_low = ceiling(quantile(n,probs = 0.025)),
#             miss_high = ceiling(quantile(n,probs = 0.975))) %>% 
#   mutate(`Setting`= "inpatient visit")

# out_miss_n <- bind_rows(out1,out2)
out_miss_n <- out1

# Compute Percentage Missed Visit Days (i.e. of all % of all missed visit days 
# during the delay window [i.e., observed - expected] days deemed to be missed 
# opportunities [here the denominator is the same across settings and since patients 
# can encounter more than one setting in a given day, these %’s will sum up to >= 100%]): 

out_miss_frac <- setting_counts1 %>% 
  mutate_at(vars(outpatient,ed,other,inpatient),~100*./n) %>% 
  select(outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean = mean(value),
            frac_median = median(value),
            frac_low = quantile(value,probs = 0.025),
            frac_high = quantile(value,probs = 0.975))


# Compute % of all missed opportunities by setting (here the %’s will sum to 100) 
out_miss_frac_miss_opps <- setting_counts1 %>% 
  select(-n) %>% 
  mutate(n = outpatient+ed+other+inpatient) %>% 
  mutate_at(vars(outpatient,ed,other,inpatient),~100*./n) %>% 
  select(outpatient,ed,other,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean_all_miss_opp = mean(value),
            frac_median_all_miss_opp = median(value),
            frac_low_all_miss_opp = quantile(value,probs = 0.025),
            frac_high_all_miss_opp = quantile(value,probs = 0.975))


### Assemble table
# Final Setting Count Table
setting_table_new_raw_all <- tibble(Setting = c("outpatient","ed","other","inpatient")) %>% 
  left_join(out_index_n, by = "Setting") %>% 
  left_join(out_pot_opps, by = "Setting") %>% 
  left_join(out_miss_n, by = "Setting") %>% 
  left_join(out_miss_frac, by = "Setting") %>% 
  left_join(out_miss_frac_miss_opps, by = "Setting")

# Output Table for Reporting
setting_counts_index_by_setting_all_new <- setting_table_new_raw_all %>% 
  mutate(`% of all Index Locations` = round(percent_index, 2)) %>% 
  mutate(`Number of Missed Opportunities` = paste0(miss_mean, " (", miss_low, "-", miss_high, ")")) %>% 
  mutate(`% of Missed opportunities` = paste0(round(frac_mean_all_miss_opp,2), " (", round(frac_low_all_miss_opp,2), "-", round(frac_high_all_miss_opp,2), ")")) %>% 
  mutate(`% of Missed Visit Days` = paste0(round(frac_mean,2), " (", round(frac_low,2), "-", round(frac_high,2), ")")) %>% 
  select(Setting,
         `Index Count`=index_n, 
         `% of all Index Locations`,
         `Potential Missed Opportunity Visit Days`=potential_opps,
         `Number of Missed Opportunities`, 
         `% of Missed opportunities`,
         `% of Missed Visit Days`) %>% 
  mutate(`% of Missed Visit Days` = ifelse(Setting == "inpatient visit", "", `% of Missed Visit Days`),
         `% of Missed opportunities` = ifelse(Setting == "inpatient visit", "", `% of Missed opportunities`))



### Save Output ----------------------------------------------------------------

save(agg_stats_all, agg_stats_ssd,
     miss_bins_all, miss_bins_ssd,
     dur_bins_all, dur_bins_ssd, 
     setting_counts_ssd, setting_counts_all,
     setting_counts_index_by_stdplac,
     setting_counts_index_by_setting,
     setting_counts_index_by_setting_ssd_new,
     setting_counts_index_by_setting_all_new,
     setting_table_new_raw_ssd,
     setting_table_new_raw_all,
     file = paste0(sim_out_path,"aggregated_sim_results.RData"))