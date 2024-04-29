## Questions ##

######################
#### Prepare Data ####
######################

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
proj_name <- "dengue"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

data_in_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

#### Load Data  ----------------------------------------------------------------
# Load in simulation results
load(paste0(sim_out_path,"sim_res_ssd.RData"))
load(paste0(sim_out_path,"sim_res_all.RData"))
load(paste0(sim_out_path,"sim_obs_reduced.RData"))
load(paste0(delay_params$base_path,"/delay_results/all_dx_visits.RData"))
load(paste0(delay_params$base_path,"/delay_results/delay_tm.RData"))
load(paste0(delay_params$out_path,"index_cases.RData"))
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

# update time map
tm <- tm %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

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
  distinct(obs,outpatient,ed,obs_stay,inpatient)

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

################################
#### Simulation_functions.R ####
################################

# generate_setting_counts_2 same as your generate_setting_counts but can handle the outer bootstrap

sim_res <- list(sim_res_obs= sim_res_ssd %>% mutate(trial = row_number()) %>% 
                  unnest(cols = c(res)),
                sim_res_miss_tab = sim_res_ssd %>% mutate(trial = row_number()) %>%
                  distinct(trial))

bootstrap_data <- boot_data %>% select(boot_trial, n_miss = n_miss_ssd)

# Step 1 take the observed time map (obs_tm) and inner join with the sim_res_ssd by obs
# this will give us the time map for the selected ssd visits for each trial. Then we do
# another inner join with tm_stdplac by obs. Note tm_stdplac can have multple observations 
# for a given obs for different stdplac for a given patient on a given day. 

tmp <- obs_tm %>%
  inner_join(sim_res$sim_res_obs, by = "obs") %>% 
  inner_join(distinct(tm_stdplac,obs,stdplac),by = join_by(obs),relationship = "many-to-many")

# Step 2 count stdplac by trial
tmp2 <- tmp %>%
  group_by(trial) %>%
  count(stdplac) %>%
  ungroup()

# Step 3 it is possible that some trials did not select a particular type of stdplac
# so we will need to add a 0 for those specific trials and stdplac so that when
# we average across trials the denominator is correct. To do this create 
# a tibble of every stdplac for each trial that we will left_join later

tmp_fill <- tmp2 %>%
  distinct(stdplac) %>%
  mutate(trial = map(stdplac,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
  unnest(trial)

# Step 4 summarize stdplac missed visits across trials

stdplac_res <- tmp2 %>%
  left_join(tmp_fill,.,by = c("stdplac", "trial")) %>%
  mutate(n=replace_na(n,0)) %>%
  group_by(trial) %>%
  mutate(pct = 100*n/sum(n)) %>%
  ungroup(trial) %>% 
  inner_join(sim_res$sim_res_obs %>% distinct(boot_trial, trial) %>% 
               inner_join(bootstrap_data, by = "boot_trial") %>% 
               distinct(trial, n_miss), by = "trial") %>% 
  mutate(pct_miss_vis_days = n/n_miss*100) %>% 
  group_by(stdplac) %>%
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
  select(stdplac, everything()) %>% 
  left_join(smallDB::stdplac_labels, by = "stdplac") 

# Step 4 now repeat for setting. Using tmp to summarise setting may be a problem
# because remember in Step 1 we since there are multiple observations per 
# obs for stdplac when we did the inner join we below up the number of observations.
# Thus, we have increased the number of outpatient, inpatient, ed, and obs_stay visits.
# If we use the following instead we wont have that problem
# tmp <- obs_tm %>%
#   inner_join(sim_res$sim_res_obs, by = "obs") 

tmp2 <- tmp %>%
  group_by(trial) %>%
  summarise(outpatient = sum(outpatient),
            ed = sum(ed),
            obs_stay = sum(obs_stay),
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

################################
#### sim_report_functions.R ####
################################

sim_res <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  unnest(cols = c(res))

tm_data = tm
bootstrap_data = boot_data %>% select(boot_trial, boot_sample)
sim_res_data = sim_res
sim_res_sim_obs_data = sim_obs_reduced %>% mutate(days_since_index = -period)
tm_stdplac_data = tm_stdplac

# generate_setting_counts_2-----------------------------------------------------

# similar to your generate_stdplac_setting_counts but can handle the outer bootstrap

# We really should get the same estimates of number of % misses as above, but we are not because we 
# are improperly average over trials by not including trials where as particular std_pla
# was not included so denominator is smaller so estimates are larger.
# Another issue is this step inner_join(tmp_index_counts,by = join_by(label)). 
# So in a give trial if there was a ssd visit at a particular stdplac but 
# no index visits at that stdplac these observations were excluded. So I fixed this in my funciton

#comment in line 150 why not exclude std_plac 81,41,42 like tou did for the index counts
#comment in line 152 why not just join by obs instead of by patient_id, svcdate, days_since_index 

# Step 1 get the number of total patients
total_patients <- filter(tm_data,days_since_index==0) %>% 
  distinct(patient_id) %>% 
  nrow()

# Step 2 get index counts by trial and stdplac 
index_stdplac_counts_temp <- tm_data %>% 
  filter(days_since_index==0) %>%
  inner_join(tm_stdplac_data,by = join_by(patient_id, svcdate, days_since_index)) %>% 
  filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
  left_join(smallDB::stdplac_labels, by = "stdplac") %>%   
  distinct(patient_id,svcdate,label) %>% 
  inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
  group_by(boot_trial) %>% 
  count(label,name = "index_visits") %>% 
  ungroup()

# Step 3 aggregate across trials 

# need to add back in 0 for trials were a particular stdplac did not occur
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

# Step 4 compute ssd counts by stdplac
# Note: I a given trial if a particular stdplac was not selected as and ssd visit
# we know nothing about its duration of delay, so this was given a value of NA
# and NA values were removed when summarizing over trials. Similarily, there 
# could be situations were for a given trial a particular stdplac was not selected as and ssd visit
# and there were no index visit with that stplac. Thus the estimate of pct_opp_missed would be 
# NaN (0/0) so these values were left as NaN and removed when summarizing over trials.

tmp  <- sim_res_data %>% 
  inner_join(sim_res_sim_obs_data, by = "obs") %>% 
  inner_join(tm_stdplac_data %>% filter(!is.na(obs)) %>% select(obs, stdplac), by = "obs",relationship = "many-to-many") %>% # na obs reflects non select sim obs
  filter(!(stdplac %in% c(81,41,42))) %>% 
  left_join(smallDB::stdplac_labels, by = "stdplac") 

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

# generate_setting_counts_3-----------------------------------------------------

# similar to your generate_main_setting_counts but can handle the outer bootstrap

# We really should get the same estimates of number of % misses as above, but we are not because here
# we are not blowing up the ed, outpatient, inpatient, obs_stays by stplac. So counts will be 
# much lower. 

# Step 1 get counts of patients

total_patients <- filter(tm_data,days_since_index==0) %>% 
  distinct(patient_id) %>% 
  nrow()

# Step 2 get index counts by trial and setting 
index_setting_counts_temp <- tm_data %>% 
  filter(days_since_index==0) %>%
  inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
  group_by(boot_trial) %>% 
  summarise(outpatient = sum(outpatient),
            ed = sum(ed),
            obs_stay= sum(obs_stay),
            inpatient = sum(inpatient)) %>% 
  ungroup() %>% 
  gather(key = Setting, value = n, -boot_trial)

# Step 3 aggregate across trials 

# need to add back in 0 for trials were a particular setting did not occur
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

# Step 4 compute ssd counts by setting
# Note: I a given trial if a particular setting was not selected as and ssd visit
# we know nothing about its duration of delay, so this was given a value of NA
# and NA values were removed when summarizing over trials. Similarily, there 
# could be situations were for a given trial a particular setting was not selected as and ssd visit
# and there were no index visit with that setting Thus the estimate of pct_opp_missed would be 
# NaN (0/0) so these values were left as NaN and removed when summarizing over trials.

tmp  <- sim_res_data %>% 
  inner_join(sim_res_sim_obs_data, by = "obs") %>% 
  inner_join(tm_data, by = join_by(patient_id, days_since_index)) %>% 
  select(trial, days_since_index, outpatient, ed, obs_stay, inpatient) %>% 
  gather(key = Setting, value = selected, -trial, -days_since_index) %>% 
  filter(selected>0)

tmp_fill <- tmp %>%
  distinct(Setting) %>%
  mutate(trial = map(Setting,~1:max(sim_res_data$trial))) %>%
  unnest(trial)

sim_res_setting_counts <- tmp %>% 
  group_by(trial,Setting) %>% 
  summarise(n=sum(selected),
            duration = mean(-days_since_index)) %>%  
  ungroup() %>% 
  left_join(tmp_fill,.,by = c("Setting", "trial")) %>%
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

sim_res_setting_counts %>% 
  full_join(index_setting_counts,by = "setting_type")