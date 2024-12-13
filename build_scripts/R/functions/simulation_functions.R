
# run single trial of simulation
run_sim_draw_visits <- function(sim_data,overall_tm){
  tmp <- run_sim_draw_simple(sim_data)
  
  
  res1 <- distinct(tmp,obs)
  
  tmp2 <- tmp %>%
    group_by(patient_id) %>%
    summarise(duration=max(period),
              n_miss=n())
  
  count_duration <- function(val){
    tmp2 %>%
      summarise(n=sum(duration>=val))
  }
  
  res2 <- tmp2  %>%
    summarise(n_pat = n(),
              mean_dur = mean(duration),
              median_dur = median(duration),
              mean_n_miss = mean(n_miss),
              median_n_miss = median(n_miss))
  
  res3 <- tibble(duration_bin = sim_data$duration_bins) %>%
    mutate(res = map(duration_bin,count_duration)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp2),
           pct_all = 100*n/sim_data$total_patients)
  
  
  count_miss <- function(val){
    tmp2 %>%
      summarise(n = sum(n_miss>=val))
  }
  
  res4 <- tibble(miss_bin = sim_data$miss_bins) %>%
    mutate(res = map(miss_bin,count_miss)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp2),
           pct_all = 100*n/sim_data$total_patients)
  
  return(list(obs = res1,
              stats = res2,
              dur_tab = res3,
              miss_tab = res4))
}

# run multiple trials of simulation
run_sim <- function(sim_data,n_trials=1000, n){
  
  cluster <- parallel::makeCluster(35)
  
  parallel::clusterCall(cluster, function() library(tidyverse))
  parallel::clusterCall(cluster, function() library(delaySim))
  parallel::clusterExport(cluster,c("run_sim_draw_visits","sim_data"),
                          envir=environment())
  
  
  sim_res <- parallel::parLapply(cl = cluster,
                                 1:n_trials,
                                 function(x){run_sim_draw_visits(sim_data)})
  
  parallel::stopCluster(cluster)
  
  
  ## Aggregate sim results
  sim_res_obs <- map2(1:n_trials,sim_res,~mutate(.y$obs,trial=.x)) %>% bind_rows()
  sim_res_stats <- map2(1:n_trials,sim_res,~mutate(.y$stats,trial=.x)) %>% bind_rows()
  sim_res_dur_tab <- map2(1:n_trials,sim_res,~mutate(.y$dur_tab,trial=.x)) %>% bind_rows()
  sim_res_miss_tab <- map2(1:n_trials,sim_res,~mutate(.y$miss_tab,trial=.x)) %>% bind_rows()
  
  return(list(sim_res_obs = sim_res_obs,
              sim_res_stats = sim_res_stats,
              sim_res_dur_tab = sim_res_dur_tab,
              sim_res_miss_tab = sim_res_miss_tab))
}

# compute simulation results to export across simulation trials
compute_export_stats <- function(sim_res,n_patients){
  ## Overall Stats ---------------------------------------------------------------
  
  sim_res_stats_means <- sim_res$sim_res_stats %>%
    summarise_all(funs(mean))
  
  sim_res_stats_upper <- sim_res$sim_res_stats %>%
    summarise_all(~quantile(.,probs = 0.975))
  
  sim_res_stats_lower <- sim_res$sim_res_stats %>%
    summarise_all(~quantile(.,probs = 0.025))
  
  sim_res_stats_means <- sim_res_stats_means %>%
    gather(key = key, value = mean) %>%
    inner_join(sim_res_stats_lower %>%
                 gather(key = key, value = low),
               by = "key") %>%
    inner_join(sim_res_stats_upper %>%
                 gather(key = key, value = high),
               by = "key") %>%
    filter(key != "trial") %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")"))
  
  # add in percent missed
  tmp <- sim_res_stats_means %>%
    filter(key == "n_pat") %>%
    mutate_at(vars(mean:high), ~100*./n_patients) %>%
    mutate(out = paste0(round(mean,2)," (",round(low,2),"-",round(high,2),")")) %>%
    mutate(key = "pct_miss")
  
  sim_res_stats_means <- bind_rows(tmp,sim_res_stats_means)
  
  rm(sim_res_stats_lower,sim_res_stats_upper,tmp)
  
  
  ## Duration stats --------------------------------------------------------------
  sim_res_dur_tab_means <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),funs(mean))
  
  
  sim_res_dur_tab_upper <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.975))
  
  sim_res_dur_tab_lower <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.025))
  
  
  sim_res_dur_tab_means <- sim_res_dur_tab_means %>%
    gather(key = key,value = mean, -duration_bin) %>%
    inner_join(sim_res_dur_tab_lower %>%
                 gather(key = key,value = low, -duration_bin),
               by = c("duration_bin", "key")) %>%
    inner_join(sim_res_dur_tab_upper %>%
                 gather(key = key,value = high, -duration_bin),
               by = c("duration_bin", "key")) %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")")) %>%
    select(duration_bin,key,out) %>%
    spread(key = key, value = out)
  
  rm(sim_res_dur_tab_lower,sim_res_dur_tab_upper)
  
  
  ## Miss stats ------------------------------------------------------------------
  sim_res_miss_tab_means <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),funs(mean))
  
  sim_res_miss_tab_upper <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.975))
  
  sim_res_miss_tab_lower <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.025))
  
  sim_res_miss_tab_means <- sim_res_miss_tab_means %>%
    gather(key = key,value = mean, -miss_bin) %>%
    inner_join(sim_res_miss_tab_lower %>%
                 gather(key = key,value = low, -miss_bin),
               by = c("miss_bin", "key")) %>%
    inner_join(sim_res_miss_tab_upper %>%
                 gather(key = key,value = high, -miss_bin),
               by = c("miss_bin", "key")) %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")")) %>%
    select(miss_bin,key,out) %>%
    spread(key = key, value = out)
  
  rm(sim_res_miss_tab_lower,sim_res_miss_tab_upper)
  
  return(list(sim_res_stats = sim_res_stats_means,
              sim_res_dur_tab = sim_res_dur_tab_means,
              sim_res_miss_tab = sim_res_miss_tab_means))
}

# generate counts of misses in each setting
generate_setting_counts <- function(obs_tm,sim_res,tm_stdplac){
  
  tmp <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>% 
    inner_join(distinct(tm_stdplac,obs,stdplac),by = join_by(obs),relationship = "many-to-many")
  
  tmp2 <- tmp %>%
    group_by(trial) %>%
    count(stdplac) %>%
    ungroup()
  
  tmp_fill <- tmp2 %>%
    distinct(stdplac) %>%
    mutate(trial = map(stdplac,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
    unnest(trial)
  
  stdplac_res <- tmp2 %>%
    left_join(tmp_fill,.,by = c("stdplac", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    group_by(stdplac) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975)) %>%
    mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
           pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    select(stdplac,n,pct,everything())
  
  
  tmp2 <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>%
    group_by(trial) %>%
    summarise(outpatient = sum(outpatient),
              ed = sum(ed),
              obs_stay = sum(obs_stay),
              inpatient = sum(inpatient)) %>% 
    ungroup() %>% 
    gather(key = setting_type, value=n, -trial)
  
  setting_type_res <- tmp2 %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    group_by(setting_type) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975)) %>%
    mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
           pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    select(setting_type,n,pct,everything())
  
  
  return(list(stdplac_res = stdplac_res,
              setting_type_res = setting_type_res))
}


# generate counts of misses in each setting version 2
generate_setting_counts_2 <- function(obs_tm, bootstrap_data, sim_res, tm_stdplac){
  
  tmp <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>% 
    inner_join(distinct(tm_stdplac,obs,stdplac),by = join_by(obs),relationship = "many-to-many")
  
  tmp2 <- tmp %>%
    group_by(trial) %>%
    count(stdplac) %>%
    ungroup()
  
  tmp_fill <- tmp2 %>%
    distinct(stdplac) %>%
    mutate(trial = map(stdplac,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
    unnest(trial)
  
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
  
  
  tmp2 <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs") %>%
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


# run interior simulation of bootstrap trial
run_boot_trials <- function(boot_ids,tm_data,miss_bins,trials,delay_params,n_patients,n_trials){
  
  tmp_sim_data <- list(time_map = boot_ids %>% 
                         inner_join(tm_data, by = "patient_id") %>% 
                         mutate(miss_ind=1L),
                       miss_bins_visits = miss_bins %>% 
                         mutate(num_miss = round(num_miss,0)),
                       change_point = delay_params$cp,
                       total_patients = n_patients,
                       miss_bins = delay_params$miss_bins,
                       duration_bins = delay_params$duration_bins)
  
  tibble(sim_trial = 1:n_trials) %>% 
    mutate(res = map(sim_trial,~run_sim_draw_simple(tmp_sim_data) %>% 
                       distinct(obs,boot_id)))
  
}