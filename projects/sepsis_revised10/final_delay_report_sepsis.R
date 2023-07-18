

library(tidyverse)
library(bit64)


cond_name <- "sepsis_revised10"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)

# final models that were run
models <- tibble(cp = delay_params$cp) %>% 
  mutate(model = map(cp,~delay_params$final_model)) %>% 
  unnest(model)

sim_out_path <- paste0(delay_params$out_path,"sim_results/")


################
#### Remove ####
################

# run markdown report



######################
#### Trends Info #####
######################

load(paste0(sim_out_path,"trends/fit_trends.RData"))
all_vis_count_fitted %>% 
  rename(tot_miss=num_miss) %>% 
  unnest(counts)


trends_out <- list()


##########################################
#### Aggregated Simulation Statistics ####
##########################################


# Loop over models




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