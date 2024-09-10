
rm(list = ls())
library(tidyverse)

cond_name <- "dengue"

# define change_point
cp <- 17

###################
#### Load Data ####
###################

load("/Shared/Statepi_Diagnosis/projects/dengue/sim_results/sim_res_ssd.RData")
# sim_res_ssd

load("/Shared/Statepi_Diagnosis/prelim_results/dengue/delay_results/sim_obs.RData")
# sim_obs

load("/Shared/Statepi_Diagnosis/prelim_results/dengue/delay_results/delay_tm.RData")
# tm

load("/Shared/Statepi_Diagnosis/prelim_results/dengue/delay_results/caseids.RData")
# caseids

load("/Shared/Statepi_Diagnosis/prelim_results/dengue/delay_results/delay_tm.RData")
# tm

load("/Shared/Statepi_Diagnosis/prelim_results/dengue/delay_results/dx_visits.RData")
# all_dx_visits

ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>% 
  mutate(dx = code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
  select(dx,dx_ver)


###############################
#### Compute Table Metrics ####
###############################

### Compute Index Setting Counts -----------------------------------------------

out_index_n <- tm %>% 
  filter(days_since_index==0) %>% 
  distinct(patient_id, outpatient,ed,inpatient,obs_stay) %>% 
  gather(key = Setting, value = value, -patient_id) %>% 
  group_by(Setting) %>% 
  summarise(index_n = sum(value))

out_index_n <- bind_rows(out_index_n,
                         filter(out_index_n,Setting=="inpatient") %>% 
                           mutate(Setting = "inpatient visit"))

### Compute Potential Opportunity Counts ---------------------------------------

# Count for visit days
out_pot_opps <- all_dx_visits %>% 
  inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  filter(between(days_since_index,-cp+1,-1)) %>% 
  distinct(obs,outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value, -obs) %>%
  group_by(Setting) %>% 
  summarise(potential_opps = sum(value))

# build counts for inpatient stays
tmp <- all_dx_visits %>% 
  inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
  inner_join(distinct(caseids,obs,caseid,patient_id),by = join_by(patient_id, obs)) %>% 
  filter(between(days_since_index,-cp+1,-1)) %>% 
  distinct(patient_id,caseid) %>% 
  count(name = "potential_opps") %>% 
  mutate(Setting = "inpatient visit")

# combine daily settings and inpatient visits
out_pot_opps <- bind_rows(out_pot_opps,tmp)

### Compute Setting Miss Counts ------------------------------------------------

obs_tm <- sim_obs %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,obs_stay,inpatient)

# # merge observations into time map to extract visit types
tmp <- sim_res_ssd %>% 
  mutate(res = map(res,~inner_join(.,obs_tm,by = "obs")))

# get setting counts for each trial
setting_counts1 <- tmp %>% 
  mutate(n = map_int(res,nrow)) %>% 
  mutate(outpatient = map_int(res,~sum(.$outpatient))) %>% 
  mutate(ed = map_int(res,~sum(.$ed))) %>% 
  mutate(obs_stay = map_int(res,~sum(.$obs_stay))) %>% 
  mutate(inpatient = map_int(res,~sum(.$inpatient))) %>% 
  select(sim_trial,boot_trial,n:inpatient)

# compute miss counts for outpatient, ed, inpatient days and obs_stay
out1 <- setting_counts1 %>% 
  select(-n) %>% 
  gather(key = Setting, value = n, -sim_trial,-boot_trial) %>% 
  group_by(Setting) %>% 
  summarise(miss_mean = round(mean(n),0),
            miss_median = round(median(n),0),
            miss_low = round(quantile(n,probs = 0.025),0),
            miss_high = round(quantile(n,probs = 0.975),0))

# compute miss counts for inpatient visits

# merge caseids into simulation results
tmp <- sim_res_ssd %>% 
  unnest(res) %>% 
  inner_join(distinct(caseids,obs,caseid,patient_id), by = "obs") %>% 
  distinct(sim_trial,boot_trial,patient_id,boot_id,caseid)

# compute counts for inpatient stays
out2 <- tmp %>% 
  count(sim_trial,boot_trial) %>% 
  summarise(miss_mean = round(mean(n),0),
            miss_median = round(median(n),0),
            miss_low = round(quantile(n,probs = 0.025),0),
            miss_high = round(quantile(n,probs = 0.975),0)) %>% 
  mutate(`Setting`= "inpatient visit")

out_miss_n <- bind_rows(out1,out2)

### Compute Percentage Miss ----------------------------------------------------

out_miss_frac <- setting_counts1 %>% 
  mutate_at(vars(outpatient,ed,obs_stay,inpatient),~round(100*./n,2)) %>% 
  select(outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean = mean(value),
            frac_median = median(value),
            frac_low = round(quantile(value,probs = 0.025),2),
            frac_high = round(quantile(value,probs = 0.975),2))


#####################################
#### Assemble Final Output Table ####
#####################################

# Final Setting Count Table
setting_table <- tibble(Setting = c("outpatient","ed","obs_stay","inpatient","inpatient visit")) %>% 
  left_join(out_index_n, by = "Setting") %>% 
  left_join(out_pot_opps, by = "Setting") %>% 
  left_join(out_miss_n, by = "Setting") %>% 
  left_join(out_miss_frac, by = "Setting")

# Output Table for Reporting
out_table <- setting_table %>% 
  mutate(`Miss Counts` = paste0(miss_mean, " (", miss_low, "-", miss_high, ")")) %>% 
  mutate(`Miss Percentage` = paste0(round(frac_mean,2), " (", frac_low, "-", frac_high, ")")) %>% 
  select(Setting,`Index Count`=index_n, `Potential Opportunities`=potential_opps,
         `Miss Counts`, `Miss Percentage`) %>% 
  mutate(`Miss Percentage` = ifelse(Setting == "inpatient visit", "", `Miss Percentage`))

