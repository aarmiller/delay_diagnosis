
# generate counts by stdplac
generate_setting_counts <- function(tm_data,sim_res_data,sim_res_sim_obs_data){
  
  # NOTE: Change this to update groupings
  new_stdplac_group <- tribble(~stdplac,~stdplac_group,
                               2, "Telehealth",
                               11, "Office / Clinic",
                               17, "Office / Clinic",
                               72, "Office / Clinic",
                               49, "Office / Clinic",
                               95, "Office / Clinic",
                               19, "Outpatient Hospital (Off Campus)",
                               20, "Urgent Care",
                               21, "Hospital (On Campus)",
                               22, "Hospital (On Campus)",
                               31, "Nursing Facility",
                               32, "Nursing Facility")
  
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()
  
  index_stdplac_counts <- tm_data %>% 
    filter(days_since_index==0) %>%
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(patient_id,stdplac_group) %>% 
    count(stdplac_group,name = "index_count") %>% 
    mutate(index_pct1 = 100*index_count/sum(index_count),   # percent of total index locations
           index_pct2 = 100*index_count/total_patients)     # percent of total patients %>% # percent of total index locations w/o Other
  
  ## get observation timemap -----------------------------------------------------
  obs_tm <- tm_data %>%
    inner_join(sim_res_sim_obs_data, by = c("patient_id", "days_since_index")) %>% 
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% 
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(obs,days_since_index,patient_id,stdplac_group) 
  
  # join obs time map into simulation results
  tmp <- sim_res_data %>% 
    inner_join(obs_tm, by = "obs")
  
  sim_res_setting_counts <- tmp %>% 
    distinct(trial,patient_id,days_since_index,stdplac_group) %>% 
    group_by(trial,stdplac_group) %>% 
    summarise(n=n(),
              duration = mean(-days_since_index)) %>%  
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    inner_join(select(index_stdplac_counts,stdplac_group,index_count),by = "stdplac_group") %>% 
    mutate(pct_opp_missed = 100*n/(n+index_count)) %>% 
    group_by(stdplac_group) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              dur_mean = mean(duration),
              dur_low = quantile(duration,probs = 0.025),
              dur_high = quantile(duration,probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              pct_opp_missed_mean = mean(pct_opp_missed),
              pct_opp_missed_low = quantile(pct_opp_missed,probs = 0.025),
              pct_opp_missed_high = quantile(pct_opp_missed,probs = 0.975)) 
  
  sim_res_setting_counts %>% 
    mutate_at(vars(n_mean:pct_opp_missed_high),~round(.,2)) %>% 
    mutate(n = paste0(n_mean," (",n_low,"-",n_high,")"),
           dur = paste0(dur_mean," (",dur_low,"-",dur_high,")"),
           pct_opp = paste0(pct_opp_mean," (",pct_opp_low,"-",pct_opp_high,")"),
           pct_opp_missed = paste0(pct_opp_missed_mean," (",pct_opp_missed_low,"-",pct_opp_missed_high,")")) %>% 
    select(stdplac_group,n,dur,pct_opp,pct_opp_missed) %>% 
    left_join(index_stdplac_counts,by = "stdplac_group")
  
}

# generate counts by stdplac version 2

generate_setting_counts_2 <- function(tm_data,sim_res_data, bootstrap_data, sim_res_sim_obs_data){
  
  # NOTE: Change this to update groupings
  new_stdplac_group <- tribble(~stdplac,~stdplac_group,
                               2, "Telehealth",
                               11, "Office / Clinic",
                               17, "Office / Clinic",
                               72, "Office / Clinic",
                               49, "Office / Clinic",
                               95, "Office / Clinic",
                               19, "Outpatient Hospital (Off Campus)",
                               20, "Urgent Care",
                               21, "Hospital (On Campus)",
                               22, "Hospital (On Campus)",
                               31, "Nursing Facility",
                               32, "Nursing Facility")
  
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()
  
  index_stdplac_counts_temp <- tm_data %>% 
    filter(days_since_index==0) %>%
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(patient_id,stdplac_group) %>% 
    inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
    group_by(boot_trial, stdplac_group) %>% 
    summarise(n= n()) %>% 
    ungroup()
  
  index_stdplac_counts <- index_stdplac_counts_temp %>% 
    group_by(boot_trial) %>% 
    mutate(index_pct1_temp = 100*n/sum(n),
           index_pct2_temp = 100*n/total_patients) %>% 
    ungroup() %>% 
    group_by(stdplac_group) %>% 
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
  obs_tm <- tm_data %>%
    inner_join(sim_res_sim_obs_data, by = c("patient_id", "days_since_index")) %>% 
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% 
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(obs,days_since_index,patient_id,stdplac_group) 
  
  # join obs time map into simulation results
  tmp <- sim_res_data %>% 
    inner_join(obs_tm, by = "obs") 
  
  sim_res_setting_counts <- tmp %>% 
    distinct(trial,patient_id,days_since_index,stdplac_group) %>% 
    group_by(trial,stdplac_group) %>% 
    summarise(n=n(),
              duration = mean(-days_since_index)) %>%  
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    inner_join(sim_res_data %>% distinct(trial, boot_trial) %>% 
                 inner_join(index_stdplac_counts_temp %>% mutate(n_index = n), by = "boot_trial") %>% 
                 distinct(trial, stdplac_group, n_index), by = c("trial", "stdplac_group")) %>% 
    mutate(pct_opp_missed = 100*n/(n+n_index)) %>% 
    group_by(stdplac_group) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              dur_mean = mean(duration),
              dur_low = quantile(duration,probs = 0.025),
              dur_high = quantile(duration,probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              pct_opp_missed_mean = mean(pct_opp_missed),
              pct_opp_missed_low = quantile(pct_opp_missed,probs = 0.025),
              pct_opp_missed_high = quantile(pct_opp_missed,probs = 0.975)) 
  
  sim_res_setting_counts %>% 
    left_join(index_stdplac_counts,by = "stdplac_group")
  
}


# generate counts by setting
generate_setting_counts_3 <- function(tm_data, sim_res_data, bootstrap_data, sim_res_sim_obs_data){
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()
  
  index_setting_counts_temp <- tm_data %>% 
    filter(days_since_index==0) %>%
    filter(setting_type!=4) %>% 
    mutate(setting_type = setting_type_labels(setting_type)) %>%
    distinct(patient_id,setting_type) %>% 
    inner_join(bootstrap_data %>% unnest(boot_sample), by = "patient_id") %>% 
    group_by(boot_trial, setting_type) %>% 
    summarise(n= n()) %>% 
    ungroup()
  
  index_setting_counts <- index_setting_counts_temp %>% 
    group_by(boot_trial) %>% 
    mutate(index_pct1_temp = 100*n/sum(n),
           index_pct2_temp = 100*n/total_patients) %>% 
    ungroup() %>% 
    group_by(setting_type) %>% 
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
  obs_tm <- tm_data %>%
    inner_join(sim_res_sim_obs_data, by = c("patient_id", "days_since_index")) %>% 
    filter(setting_type!=4) %>% 
    mutate(setting_type = setting_type_labels(setting_type)) %>%
    distinct(obs,days_since_index,patient_id,setting_type) 
  
  # join obs time map into simulation results
  tmp <- sim_res_data %>% 
    inner_join(obs_tm, by = "obs")
  
  sim_res_setting_counts <- tmp %>% 
    distinct(trial,patient_id,days_since_index,setting_type) %>% 
    group_by(trial,setting_type) %>% 
    summarise(n=n(),
              duration = mean(-days_since_index)) %>%  
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    inner_join(sim_res_data %>% distinct(trial, boot_trial) %>% 
                 inner_join(index_setting_counts_temp %>% mutate(n_index = n), by = "boot_trial") %>% 
                 distinct(trial, setting_type, n_index), by = c("trial", "setting_type")) %>% 
    mutate(pct_opp_missed = 100*n/(n+n_index)) %>% 
    group_by(setting_type) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              dur_mean = mean(duration),
              dur_low = quantile(duration,probs = 0.025),
              dur_high = quantile(duration,probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              pct_opp_missed_mean = mean(pct_opp_missed),
              pct_opp_missed_low = quantile(pct_opp_missed,probs = 0.025),
              pct_opp_missed_high = quantile(pct_opp_missed,probs = 0.975)) 
  
  sim_res_setting_counts %>% 
    left_join(index_setting_counts,by = "setting_type")
  
}

