
rm(list = ls())

library(tidyverse)
library(bit64)

cond_name <- "measles"
load("/Volumes/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

delay_params


load(paste0("/Volumes/Statepi_Diagnosis/projects/measles/", "index_cases.RData"))

fams <- index_cases %>% 
  mutate(efamid = enrolid %/% 100) %>% 
  count(efamid) %>% 
  filter(n>1) %>% 
  distinct(efamid)

tmp <- fams %>% 
  inner_join(index_cases %>% 
               mutate(efamid = enrolid %/% 100)) %>% 
  group_by(efamid) %>% 
  mutate(fam_index = min(index_date)) %>% 
  ungroup() %>% 
  mutate(fam_first = ifelse(index_date==fam_index,1L,0L))

study_ids <- bind_rows(tmp %>% 
  filter(fam_first == 0) %>% 
  mutate(days_since_fam = index_date-fam_index) %>% 
  filter(days_since_fam<=60) %>% 
  distinct(patient_id,fam_first),
  tmp %>% 
    filter(fam_first == 1) %>% 
    distinct(patient_id,fam_first))

load("/Volumes/Statepi_Diagnosis/projects/measles/delay_window_1_14/sim_results/sim_obs_reduced.RData")



load("/Volumes/Statepi_Diagnosis/projects/measles/delay_window_1_14/sim_results/boot_data.RData")

trial_cases <- select(boot_data,boot_trial,boot_sample) %>% 
  mutate(boot_sample = map(boot_sample,~inner_join(.,study_ids,by = join_by(patient_id))))

# Check if any secondary cases had a miss

tmp1 <- sim_obs_reduced %>% 
  inner_join(study_ids) %>% 
  filter(fam_first == 0)

sim_res_ssd %>% 
  select(sim_trial,res) %>% 
  unnest(res) %>% 
  inner_join(tmp1)

# NOPE


tmp_res <- sim_res_ssd %>% 
  arrange(boot_trial,sim_trial) %>% 
  mutate(trial = row_number()) %>% 
  select(trial,everything())

# filter to the fam enrolids in a given trial
tmp_res <- tmp_res %>% 
  inner_join(trial_cases,by = join_by(boot_trial)) %>% 
  mutate(res = map2(res,boot_sample,~inner_join(.x,.y,by = join_by(boot_id))))

tmp_res %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  count(fam_first)

tmp_res <- tmp_res %>% 
  mutate(res = map(res,~inner_join(.,sim_obs_reduced, by = join_by(obs, patient_id))))

tmp_res <- tmp_res %>% 
  mutate(res = map(res,~group_by(.,boot_id) %>% 
                     summarise(delay = max(period))))

tmp_res <- tmp_res %>% 
  mutate(boot_sample = map2(boot_sample,res,~left_join(.x,.y, by = "boot_id")))


tmp <- tmp_res %>% 
  select(trial,boot_sample) %>% 
  unnest(boot_sample)

tmp <- tmp %>% 
  mutate(miss_patient = ifelse(is.na(delay),0L,1L)) %>% 
  group_by(trial,fam_first) %>% 
  summarise(n_cases = n(),
            n_misses = sum(miss_patient),
            mean_duration = mean(delay,na.rm = T),
            median_duration = median(delay,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(mean_duration = replace_na(mean_duration,0)) %>% 
  mutate(median_duration = replace_na(median_duration,0))


tmp %>% 
  group_by(fam_first) %>% 
  summarise(n_miss = mean(n_misses))

sim_obs_reduced %>% 
  inner_join(study_ids) %>% 
  count(fam_first)

tmp %>% 
  filter(fam_first == 0) %>% 
  filter(n_misses>0)
  

select(tmp_res,trial,res) %>% 
  unnest(res) %>% 
  inner_join(sim_obs_reduced)
