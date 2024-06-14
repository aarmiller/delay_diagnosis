
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
    # filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
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
    # filter(!(stdplac %in% c(81,41,42))) %>% 
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


# NEW - generate the main setting counts
generate_main_setting_counts <- function(tm_data,sim_res_data,sim_res_sim_obs_data){
  
  tmp_index_counts <- tm_data %>% 
    filter(days_since_index==0) %>%
    summarise(Outpatient = sum(outpatient),
              ED = sum(ed),
              `Observational Stay` = sum(obs_stay),
              Inpatient = sum(inpatient)) %>% 
    gather(key = Setting,value = index_count) %>% 
    mutate(index_pct1 = 100*index_count/sum(index_count),   # percent of total index locations
           index_pct2 = 100*index_count/total_patients)     # percent of total patients
  
  out <- sim_res_data %>% 
    inner_join(sim_res_sim_obs_data, by = "obs") %>% 
    inner_join(tm_data, by = join_by(patient_id, days_since_index)) %>% 
    group_by(trial) %>% 
    summarise(Outpatient = sum(outpatient),
              ED = sum(ed),
              `Observational Stay` = sum(obs_stay),
              Inpatient = sum(inpatient)) %>% 
    gather(key = Setting,value = n, -trial) %>% 
    inner_join(select(tmp_index_counts,Setting,index_count),by = join_by(Setting)) %>% 
    mutate(total_opportunities = n+index_count,
           pct_opp_miss = 100*n/total_opportunities) %>% 
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    group_by(Setting) %>%
    summarise(miss_mean = mean(n),
              miss_low = quantile(n,probs = 0.025),
              miss_high = quantile(n,probs = 0.975),
              tot_opp_mean = mean(total_opportunities),
              tot_opp_low = quantile(total_opportunities,probs = 0.025),
              tot_opp_high = quantile(total_opportunities,probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              pct_opp_miss_mean = mean(pct_opp_miss),
              pct_opp_miss_low = quantile(pct_opp_miss,probs = 0.025),
              pct_opp_miss_high = quantile(pct_opp_miss,probs = 0.975))
  
  out <- inner_join(tmp_index_counts,out, by = "Setting")
  
  return(out)
  
}

# NEW - generate the stdplac counts
generate_stdplac_setting_counts <- function(tm_data,sim_res_data,sim_res_sim_obs_data,tm_stdplac_data){
  
  tmp_index_counts <- tm_data %>% 
    filter(days_since_index==0) %>%
    inner_join(tm_stdplac_data,by = join_by(patient_id, svcdate, days_since_index)) %>% 
    # filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(smallDB::stdplac_labels, by = "stdplac") %>% 
    distinct(patient_id,svcdate,label) %>% 
    count(label,name = "index_visits") %>% 
    arrange(desc(index_visits)) 
  
  out <- sim_res_data %>% 
    inner_join(sim_res_sim_obs_data, by = "obs") %>% 
    inner_join(tm_stdplac_data,by = join_by(patient_id, days_since_index),relationship = "many-to-many") %>% 
    left_join(smallDB::stdplac_labels, by = "stdplac") %>% 
    group_by(trial) %>% 
    count(label) %>% 
    inner_join(tmp_index_counts,by = join_by(label)) %>% 
    mutate(total_opportunities = n+index_visits,
           pct_opp_miss = 100*n/total_opportunities) %>% 
    group_by(trial) %>% 
    mutate(pct_opp = 100*n/sum(n)) %>% 
    ungroup() %>% 
    group_by(label) %>%
    summarise(miss_mean = mean(n),
              miss_low = quantile(n,probs = 0.025),
              miss_high = quantile(n,probs = 0.975),
              tot_opp_mean = mean(total_opportunities),
              tot_opp_low = quantile(total_opportunities,probs = 0.025),
              tot_opp_high = quantile(total_opportunities,probs = 0.975),
              pct_opp_mean = mean(pct_opp),
              pct_opp_low = quantile(pct_opp,probs = 0.025),
              pct_opp_high = quantile(pct_opp,probs = 0.975),
              pct_opp_miss_mean = mean(pct_opp_miss),
              pct_opp_miss_low = quantile(pct_opp_miss,probs = 0.025),
              pct_opp_miss_high = quantile(pct_opp_miss,probs = 0.975)) 
  
  out <- inner_join(tmp_index_counts,out, by = "label") %>%  
    rename(`Place of Care`=label)
  
  return(out)
}

# generate counts by stdplac version 2
generate_setting_counts_2 <- function(tm_data,sim_res_data, bootstrap_data, sim_res_sim_obs_data, tm_stdplac_data){
  
  # count patients
  total_patients <- filter(tm_data,days_since_index==0) %>% 
    distinct(patient_id) %>% 
    nrow()

  index_stdplac_counts_temp <- tm_data %>% 
    filter(days_since_index==0) %>%
    inner_join(tm_stdplac_data,by = join_by(patient_id, svcdate, days_since_index)) %>% 
    # filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(smallDB::stdplac_labels, by = "stdplac") %>%   
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
    inner_join(tm_stdplac_data %>% filter(!is.na(obs)) %>% select(obs, stdplac), by = "obs",relationship = "many-to-many") %>% # na obs reflects non select sim obs
    # filter(!(stdplac %in% c(81,41,42))) %>% 
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
  
}


# generate counts by setting
generate_setting_counts_3 <- function(tm_data, sim_res_data, bootstrap_data, sim_res_sim_obs_data, tm_stdplac_data){
  
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
              obs_stay= sum(obs_stay),
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
    select(trial, days_since_index, outpatient, ed, obs_stay, inpatient)
  
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


