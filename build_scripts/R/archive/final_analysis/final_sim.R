
rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)
library(codeBuildr)
# devtools::install_github("aarmiller/codeBuildr")

setwd("/Shared/Statepi_Diagnosis/atlan")

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

# source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
# source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

args = commandArgs(trailingOnly=TRUE)

# name of condition
proj_name <- args[1]
# proj_name <- "sarcoid_skin"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

data_in_path <- paste0(delay_params$base_path,"delay_results/")
# sim_in_path <- paste0("/Volumes/AML/small_dbs/",cond_name,"/truven/enroll_restrict_365/","delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

#### Load Data  ----------------------------------------------------------------

load(paste0(data_in_path,"all_dx_visits.RData"))
load(paste0(data_in_path,"delay_tm.RData"))
load(paste0(delay_params$out_path,"index_cases.RData"))

# load alls
project_test <- codeBuildr::avail_ssd_codes() %>% 
  filter(name == proj_name) %>% nrow()

if(project_test>0){

  ssd_codes <- codeBuildr::load_ssd_codes(proj_name) %>% 
    filter(type %in% c("icd9","icd10")) %>% 
    mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
    select(dx = code,dx_ver)
  
} else {

  ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>% 
    filter(type %in% c("icd9","icd10")) %>% 
    mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
    select(dx = code,dx_ver)
}

# extract patient ids and number of patients
patient_ids <- index_cases %>% distinct(patient_id)

n_patients <- nrow(patient_ids)


#####################################
#### Compute Overall Trends/Fits ####
#####################################

# Here we compute the overall trends and number of misses outside of the bootstrap

### Compute visit counts -------------------------------------------------------

# all visits
sim_tm_all <- all_dx_visits %>%
  inner_join(patient_ids, by = "patient_id") %>% 
  mutate(period = -days_since_index) %>%
  distinct(patient_id,period,days_since_index) %>%
  inner_join(sim_obs,by = c("patient_id", "days_since_index")) 

all_vis_count <- sim_tm_all %>% 
  count(period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

# ssd visits
sim_tm_ssd <- all_dx_visits %>%
  inner_join(patient_ids, by = "patient_id") %>% 
  mutate(period = -days_since_index) %>%
  inner_join(ssd_codes,by = c("dx", "dx_ver")) %>% 
  distinct(patient_id,period,days_since_index) %>%
  inner_join(sim_obs,by = c("patient_id", "days_since_index")) 

ssd_vis_count <- sim_tm_ssd %>% 
  count(.,period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))


#### fit trends ----------------------------------------------------------------

# all visits
all_vis_count <- return_fits(data = all_vis_count,
                             model = delay_params$final_model,
                             cp = delay_params$cp,
                             periodicity = delay_params$periodicity)

# ssd visits
ssd_vis_count <- return_fits(ssd_vis_count,
                             model = delay_params$final_model,
                             cp = delay_params$cp,
                             periodicity = delay_params$periodicity)

# miss counts
# all_vis_count %>% summarise(n_miss = sum(num_miss,na.rm = T))
# ssd_vis_count %>% summarise(n_miss = sum(num_miss,na.rm = T))
save(all_vis_count, ssd_vis_count,
     file = paste0(sim_out_path,"fit_trends.RData"))

#### Export Plot ---------------------------------------------------------------

bind_rows(mutate(ssd_vis_count,group = "SSD Visits"),
          mutate(all_vis_count,group = "ALL Visits")) %>% 
  ggplot(aes(period,n)) +
  geom_line() +
  geom_line(aes(y = pred), color = "red") +
  facet_wrap(~group) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp)) +
  theme_bw()
ggsave(paste0(sim_out_path,"visit_trends.jpeg"), width = 10, height = 6)



################################
#### Prepare Bootstrap Data ####
################################

#### Setup outer bootstrap data ------------------------------------------------

# sample overall set of patients for each bootstrap
tmp <- tibble(boot_trial=1:delay_params$boot_trials) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(patient_id = sample(patient_ids$patient_id, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# tmp$boot_sample[[1]]

#### compute counts for all visits ---------------------------------------------
# For each set of patients selected pull the time map visit counts
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "patient_id") %>%
                        mutate(period = -days_since_index) %>%
                        distinct(patient_id,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("patient_id", "days_since_index"))))

tmp <- tmp %>% 
  mutate(all_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7))))

#### compute counts for ssd visits ---------------------------------------------
# for each set of patients pull the ssd counts
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "patient_id") %>%
                        mutate(period = -days_since_index) %>%
                        inner_join(ssd_codes,by = c("dx", "dx_ver")) %>%
                        distinct(patient_id,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("patient_id", "days_since_index"))))

tmp <- tmp %>% 
  mutate(ssd_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7))))

# remove the sim_tm data (no longer needed)
tmp <- tmp %>% 
  select(boot_trial,boot_sample,all_vis_count,ssd_vis_count)


### fit models -----------------------------------------------------------------

# fit trends for both all visits and ssd visits
tmp <- tmp %>% 
  mutate(all_vis_count = map(all_vis_count,~return_fits(.,
                                                        model = delay_params$final_model,
                                                        cp = delay_params$cp,
                                                        periodicity = delay_params$periodicity)),
         ssd_vis_count = map(ssd_vis_count,~return_fits(.,
                                                        model = delay_params$final_model,
                                                        cp = delay_params$cp,
                                                        periodicity = delay_params$periodicity)))
# compute number missed
tmp <- tmp %>% 
  mutate(n_miss_all = map(all_vis_count,
                          ~summarise(.,n_miss_all = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_all) %>% 
  mutate(n_miss_ssd = map(ssd_vis_count,
                          ~summarise(.,n_miss_ssd = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_ssd)

#### finalize and save data for bootstrapping ----------------------------------
boot_data <- tmp
rm(tmp)

# save bootstrap data
save(boot_data,file = paste0(sim_out_path,"boot_data.RData"))

# Compute number missed (note: need to export info in a report in future)
# boot_data %>% 
#   summarise(mean_n_miss_all = mean(n_miss_all),
#             median_n_miss_all = median(n_miss_all),
#             mean_n_miss_ssd = mean(n_miss_ssd),
#             median_n_miss_ssd = median(n_miss_ssd))



###########################################
#### Run bootstrap simulation analysis ####
###########################################

cluster <- parallel::makeCluster(30)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterCall(cluster, function() library(delaySim))
parallel::clusterExport(cluster,c("run_boot_trials","boot_data","sim_tm_all","sim_tm_ssd","delay_params","n_patients"),
                        envir=environment())

set.seed(5678)
sim_res_all <- parallel::parLapply(cl = cluster,
                                   1:delay_params$boot_trials,
                                   function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                               tm_data = sim_tm_all,
                                                               miss_bins = select(boot_data$all_vis_count[[x]],period,num_miss),
                                                               delay_params = delay_params,
                                                               n_patients = n_patients,
                                                               n_trials = delay_params$sim_trials)})

sim_res_all <- map2(sim_res_all,1:delay_params$boot_trials,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

set.seed(5678)
sim_res_ssd <- parallel::parLapply(cl = cluster,
                                   1:delay_params$boot_trials,
                                   function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                               tm_data = sim_tm_ssd,
                                                               miss_bins = select(boot_data$ssd_vis_count[[x]],period,num_miss),
                                                               delay_params = delay_params,
                                                               n_patients = n_patients,
                                                               n_trials = delay_params$sim_trials)})

sim_res_ssd <- map2(sim_res_ssd,1:delay_params$boot_trials,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

parallel::stopCluster(cluster)

rm(cluster)

save(sim_res_all,file = paste0(sim_out_path,"sim_res_all.RData"))
save(sim_res_ssd,file = paste0(sim_out_path,"sim_res_ssd.RData"))

#### Extract sim result observations -------------------------------------------

# now extract the specific sim obs values for future steps
tmp1 <- sim_res_ssd %>% unnest(res) %>% distinct(obs)
tmp2 <- sim_res_all %>% unnest(res) %>% distinct(obs)

sim_obs_reduced <- sim_obs %>% 
  inner_join(distinct(bind_rows(tmp1,tmp2)), by = "obs") %>% 
  mutate(period = -days_since_index) %>% 
  select(obs,patient_id,period)

save(sim_obs_reduced, file = paste0(sim_out_path,"sim_obs_reduced.RData"))

# cleanup
rm(list = ls()[!(ls() %in% c("delay_params","sim_out_path","sim_in_path","cond_name","n_patients"))])
gc()



# I moved the following to a separate script that we can run separately if we need to 
# change parameters for summary stats, without need to rerun simulations.

# ###################################
# #### Analyze Bootstrap Results ####
# ###################################
# 
# # Load in sumulation results
# load(paste0(sim_out_path,"sim_res_ssd.RData"))
# load(paste0(sim_out_path,"sim_res_all.RData"))
# load(paste0(sim_out_path,"sim_obs_reduced.RData"))
# 
# compute_boot_stats <- function(sim_data,sim_obs_data,delay_params,n_patients){
#   
#   tmp1 <- sim_data %>% 
#     inner_join(sim_obs_data, by = "obs") %>% 
#     group_by(patient_id,boot_id) %>%
#     summarise(duration=max(period),
#               n_miss=n())
#   
#   count_duration <- function(val){
#     tmp1 %>%
#       ungroup() %>% 
#       summarise(n=sum(duration>=val))
#   }
#   
#   res1 <- tmp1  %>%
#     ungroup() %>% 
#     summarise(n_pat = n(),
#               mean_dur = mean(duration),
#               median_dur = median(duration),
#               mean_n_miss = mean(n_miss),
#               median_n_miss = median(n_miss))
#   
#   res2 <- tibble(duration_bin = delay_params$duration_bins) %>%
#     mutate(res = map(duration_bin,count_duration)) %>%
#     unnest(res) %>%
#     mutate(pct_miss = 100*n/nrow(tmp1),
#            pct_all = 100*n/n_patients)
#   
#   count_miss <- function(val){
#     tmp1 %>%
#       ungroup() %>% 
#       summarise(n = sum(n_miss>=val))
#   }
#   
#   res3 <- tibble(miss_bin = delay_params$miss_bins) %>%
#     mutate(res = map(miss_bin,count_miss)) %>%
#     unnest(res) %>%
#     mutate(pct_miss = 100*n/nrow(tmp1),
#            pct_all = 100*n/n_patients)
#   
#   return(list(stats = res1,
#               dur_tab = res2,
#               miss_tab = res3))
#   
# }
# 
# # setup cluster 
# plan(multisession, workers = 30)
# 
# tmp1 <- sim_res_ssd %>% 
#   # slice(1:30) %>% 
#   mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
#                                                          sim_obs_reduced,
#                                                          delay_params = delay_params, 
#                                                          n_patients = n_patients)))
# 
# tmp1 <- tmp1 %>% 
#   select(sim_trial,boot_trial,sim_stats) %>% 
#   mutate(stats = map(sim_stats,~.$stats)) %>% 
#   mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>% 
#   mutate(dur_tab = map(sim_stats,~.$dur_tab))
# 
# sim_stats_ssd <- tmp1 %>% 
#   select(-sim_stats)
# 
# tmp2 <- sim_res_all %>% 
#   # slice(1:30) %>% 
#   mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
#                                                          sim_obs_reduced,
#                                                          delay_params = delay_params, 
#                                                          n_patients = n_patients)))
# 
# tmp2 <- tmp2 %>% 
#   select(sim_trial,boot_trial,sim_stats) %>% 
#   mutate(stats = map(sim_stats,~.$stats)) %>% 
#   mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>% 
#   mutate(dur_tab = map(sim_stats,~.$dur_tab))
# 
# sim_stats_all <- tmp2 %>% 
#   select(-sim_stats)
# 
# save(sim_stats_ssd,sim_stats_all,file = paste0(sim_out_path,"sim_stats.RData"))
# 
# 
# ###################################
# #### Export Simulation Results ####
# ###################################
# 
# # Note: We may want to move parts of this section to a separate markdown report script
# 
# ### Generate aggregate statistics ----------------------------------------------
# 
# aggregate_sim_stats <- function(sim_stats_data){
#   tmp1 <- sim_stats_data %>% 
#     select(sim_trial,boot_trial,stats) %>% 
#     unnest(stats) %>% 
#     summarise_at(vars(n_pat:median_n_miss),~round(mean(.),2))
#   
#   tmp2 <- sim_stats_data %>% 
#     select(sim_trial,boot_trial,stats) %>% 
#     unnest(stats) %>% 
#     summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.025)),2))
#   
#   tmp3 <- sim_stats_data %>% 
#     select(sim_trial,boot_trial,stats) %>% 
#     unnest(stats) %>% 
#     summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.975)),2))
#   
#   
#   res1 <- inner_join(gather(tmp1, key=measure, value = mean),
#                      gather(tmp2, key=measure, value = low)) %>% 
#     inner_join(gather(tmp3, key=measure, value = high)) %>% 
#     mutate(measure_out=paste0(mean," (",low,"-",high,")"))
#   
#   return(list(main_stats = res1))
# }
# 
# # generate aggregated results of:
# #   n_pat - number of patient missed
# #   mean_dur - mean duration of misses
# #   median_dur - median duration of misses
# #   mean_n_miss - mean number of misses per patient
# #   median_n_miss - median number of misses per patient
# agg_stats_ssd <- aggregate_sim_stats(sim_stats_ssd) 
# agg_stats_all <- aggregate_sim_stats(sim_stats_all)
# 
# ### Generate bins for number of misses per patient -----------------------------
# 
# # SSD Visits
# tmp1 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
#   gather(key = key, value = mean, -miss_bin) %>% 
#   mutate(mean = round(mean,2))
# 
# tmp2 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
#   gather(key = key, value = low, -miss_bin) %>% 
#   mutate(low = round(low,2))
# 
# tmp3 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
#   gather(key = key, value = high, -miss_bin) %>% 
#   mutate(high = round(high,2))
# 
# miss_bins_ssd <- inner_join(tmp1,tmp2) %>% 
#   inner_join(tmp3) %>% 
#   mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
#   select(miss_bin,key,out) %>% 
#   spread(key = key, value = out)
# 
# # All visits
# tmp1 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
#   gather(key = key, value = mean, -miss_bin) %>% 
#   mutate(mean = round(mean,2))
# 
# tmp2 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
#   gather(key = key, value = low, -miss_bin) %>% 
#   mutate(low = round(low,2))
# 
# tmp3 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,miss_tab) %>% 
#   unnest(miss_tab) %>% 
#   group_by(miss_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
#   gather(key = key, value = high, -miss_bin) %>% 
#   mutate(high = round(high,2))
# 
# miss_bins_all <- inner_join(tmp1,tmp2) %>% 
#   inner_join(tmp3) %>% 
#   mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
#   select(miss_bin,key,out) %>% 
#   spread(key = key, value = out)
# 
# 
# ### Generate bins for duration of misses ---------------------------------------
# 
# # SSD Visits
# tmp1 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
#   gather(key = key, value = mean, -duration_bin) %>% 
#   mutate(mean = round(mean,2))
# 
# tmp2 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
#   gather(key = key, value = low, -duration_bin) %>% 
#   mutate(low = round(low,2))
# 
# tmp3 <- sim_stats_ssd %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
#   gather(key = key, value = high, -duration_bin) %>% 
#   mutate(high = round(high,2))
# 
# dur_bins_ssd <- inner_join(tmp1,tmp2) %>% 
#   inner_join(tmp3) %>% 
#   mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
#   select(duration_bin,key,out) %>% 
#   spread(key = key, value = out)
# 
# # ALL Visits
# tmp1 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(mean(.),2)) %>% 
#   gather(key = key, value = mean, -duration_bin) %>% 
#   mutate(mean = round(mean,2))
# 
# tmp2 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>% 
#   gather(key = key, value = low, -duration_bin) %>% 
#   mutate(low = round(low,2))
# 
# tmp3 <- sim_stats_all %>% 
#   select(sim_trial,boot_trial,dur_tab) %>% 
#   unnest(dur_tab) %>% 
#   group_by(duration_bin) %>% 
#   summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>% 
#   gather(key = key, value = high, -duration_bin) %>% 
#   mutate(high = round(high,2))
# 
# dur_bins_all <- inner_join(tmp1,tmp2) %>% 
#   inner_join(tmp3) %>% 
#   mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
#   select(duration_bin,key,out) %>% 
#   spread(key = key, value = out)
# 
# ### Save Output ----------------------------------------------------------------
# 
# save(agg_stats_all, agg_stats_ssd,
#      miss_bins_all, miss_bins_ssd,
#      dur_bins_all, dur_bins_ssd, 
#      file = paste0(sim_out_path,"aggregated_sim_results.RData"))