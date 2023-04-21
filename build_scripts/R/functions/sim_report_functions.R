
# generate counts by setting
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
    distinct(enrolid) %>% 
    nrow()
  
  index_stdplac_counts <- tm_data %>% 
    filter(days_since_index==0) %>%
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% #remove lab (81), ambulance land (41), ambulance air/water (42)
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(enrolid,stdplac_group) %>% 
    count(stdplac_group,name = "index_count") %>% 
    mutate(index_pct1 = 100*index_count/sum(index_count),   # percent of total index locations
           index_pct2 = 100*index_count/total_patients)     # percent of total patients %>% # percent of total index locations w/o Other
  
  ## get observation timemap -----------------------------------------------------
  obs_tm <- tm_data %>%
    inner_join(sim_res_sim_obs_data, by = c("enrolid", "days_since_index")) %>% 
    filter(setting_type!=4) %>% 
    filter(!(stdplac %in% c(81,41,42))) %>% 
    left_join(new_stdplac_group, by = "stdplac") %>% 
    mutate(stdplac_group = ifelse(setting_type==2,"ED",
                                  ifelse(setting_type==5,"Inpatient",stdplac_group))) %>% 
    mutate(stdplac_group=replace_na(stdplac_group,"Other Outpatient")) %>% 
    distinct(obs,days_since_index,enrolid,stdplac_group) 
  
  # join obs time map into simulation results
  tmp <- sim_res_data %>% 
    inner_join(obs_tm, by = "obs")
  
  sim_res_setting_counts <- tmp %>% 
    distinct(trial,enrolid,days_since_index,stdplac_group) %>% 
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
