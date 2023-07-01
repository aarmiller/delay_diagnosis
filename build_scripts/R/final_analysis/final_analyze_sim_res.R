
rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)

setwd("/Shared/Statepi_Diagnosis/atlan")

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

args = commandArgs(trailingOnly=TRUE)

# name of condition
proj_name <- args[1]
# proj_name <- "sarcoid_lung"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

data_in_path <- paste0(delay_params$base_path,"delay_results/")
# sim_in_path <- paste0("/Volumes/AML/small_dbs/",proj_name,"/truven/enroll_restrict_365/","delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

######################
#### Prepare Data ####
######################

#### Load Data  ----------------------------------------------------------------
# Load in sumulation results
load(paste0(sim_out_path,"sim_res_ssd.RData"))
load(paste0(sim_out_path,"sim_res_all.RData"))
load(paste0(sim_out_path,"sim_obs_reduced.RData"))
load(paste0(delay_params$base_path,"/delay_results/all_dx_visits.RData"))
load(paste0(delay_params$base_path,"/delay_results/delay_tm.RData"))
load(paste0(delay_params$out_path,"index_cases.RData"))

# Subset to project specific patient ids
index_dx_dates <- index_dx_dates %>% inner_join(index_cases)
patient_ids <- index_cases %>% distinct(patient_id)
all_dx_visits <- all_dx_visits %>% inner_join(patient_ids)
total_patients <- nrow(patient_ids)
tm <- tm %>% inner_join(patient_ids)

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
  distinct(obs,stdplac,setting_type) %>%
  filter(setting_type!=4)

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
plan(multisession, workers = 30)

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

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")

sim_res <- list(sim_res_obs= sim_res_ssd %>% mutate(trial = row_number()) %>% 
                  unnest(cols = c(res)),
                sim_res_miss_tab = sim_res_ssd %>% mutate(trial = row_number()) %>%
                  distinct(trial))

setting_counts_ssd <- generate_setting_counts(obs_tm = obs_tm,
                                          sim_res =sim_res)

sim_res <- list(sim_res_obs= sim_res_all %>% mutate(trial = row_number()) %>% 
                  unnest(cols = c(res)),
                sim_res_miss_tab = sim_res_all %>% mutate(trial = row_number()) %>%
                  distinct(trial))
setting_counts_all <- generate_setting_counts(obs_tm = obs_tm,
                                          sim_res = sim_res)

### Setting summary index ------------------------------------------------------

source("github/delay_diagnosis/build_scripts/R/functions/sim_report_functions.R")
sim_res <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  unnest(cols = c(res))

# Generate the index counts
setting_counts_index <- generate_setting_counts(tm_data = tm,
                               sim_res_data = sim_res,
                               sim_res_sim_obs_data = sim_obs_reduced %>% mutate(days_since_index = -period))

### Rurality -------------------------------------------------------------------

load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% inner_join(patient_ids)
demo2 <- demo2 %>% inner_join(patient_ids)

rural_ids <- rural_visits %>% inner_join(patient_ids) %>% distinct(patient_id)

msa_ids <- demo2 %>% 
  filter(index_date>=dtstart & index_date<=dtend) %>% 
  filter(!is.na(msa) & msa!="0" & msa != "") %>% 
  distinct(patient_id)

no_msa_ids <- demo2 %>% 
  filter(index_date>=dtstart & index_date<=dtend) %>% 
  filter(msa=="0" | msa=="") %>%
  distinct(patient_id)

rural_no_msa_ids <- rural_ids %>% anti_join(msa_ids)

get_id_res <- function(id_set){
  sim_obs_reduced %>% mutate(days_since_index = -period) %>% 
    inner_join(id_set) %>% 
    inner_join(sim_res) %>% 
    group_by(trial,patient_id) %>% 
    summarise(duration = -min(days_since_index),
              n_miss = n()) %>% 
    group_by(trial) %>% 
    summarise(n_pat = n(),
              n_miss = mean(n_miss),
              duration = mean(duration)) %>% 
    summarise(n_pat_mean = mean(n_pat),
              n_pat_low = quantile(n_pat,probs = c(0.025)),
              n_pat_high = quantile(n_pat,probs = c(0.975)),
              n_miss_mean = mean(n_miss),
              n_miss_low = quantile(n_miss,probs = c(0.025)),
              n_miss_high = quantile(n_miss,probs = c(0.975)),
              duration_mean = mean(duration),
              duration_low = quantile(duration,probs = c(0.025)),
              duration_high = quantile(duration,probs = c(0.975))) %>% 
    mutate(pct_miss_mean = 100*n_pat_mean/nrow(id_set),
           pct_miss_low = 100*n_pat_low/nrow(id_set),
           pct_miss_high = 100*n_pat_high/nrow(id_set)) %>% 
    mutate_all(~round(.,2)) %>% 
    mutate(n_pat = paste0(n_pat_mean," (",n_pat_low,"-",n_pat_high,")"),
           pct_miss = paste0(pct_miss_mean," (",pct_miss_low,"-",pct_miss_high,")"),
           duration = paste0(duration_mean," (",duration_low,"-",duration_high,")"),
           n_miss = paste0(n_miss_mean," (",n_miss_low,"-",n_miss_high,")")) %>% 
    mutate(total_patients = nrow(id_set)) %>%
    select(total_patients,n_pat:n_miss) %>% 
    gather(key = key, value = value)
}

rural_res <- get_id_res(rural_ids) %>% 
  rename(Rural = value)

msa_res <- get_id_res(msa_ids) %>% 
  rename(`In MSA` = value)

rural_no_msa_res <- get_id_res(rural_no_msa_ids) %>% 
  rename(`Rural No MSA` = value)

no_msa_res <-  get_id_res(no_msa_ids) %>% 
  rename(`No MSA` = value)

location_counts <- msa_res %>% 
  inner_join(no_msa_res) %>% 
  inner_join(rural_res) %>% 
  inner_join(rural_no_msa_res) %>% 
  rename(Measure = key)

### Save Output ----------------------------------------------------------------

save(agg_stats_all, agg_stats_ssd,
     miss_bins_all, miss_bins_ssd,
     dur_bins_all, dur_bins_ssd, 
     setting_counts_ssd, setting_counts_all,
     setting_counts_index,
     location_counts,
     file = paste0(sim_out_path,"aggregated_sim_results.RData"))