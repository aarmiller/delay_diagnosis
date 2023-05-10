


rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "tb"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

sim_in_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")
# sim_in_path <- paste0("/Volumes/AML/small_dbs/",cond_name,"/truven/enroll_restrict_365/","delay_results/")
sim_out_path <- paste0(delay_params$base_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

#### Load Data  ----------------------------------------------------------------

load(paste0(sim_in_path,"all_dx_visits.RData"))
load(paste0(sim_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver)

# extract patient ids and number of patients
patient_ids <- tm %>% 
  distinct(enrolid)

n_patients <- nrow(patient_ids)


########################################
###### Apply Filtering Criteria ########
########################################

# Apply any criteria to restrict to final study population here


#####################################
#### Compute Overall Trends/Fits ####
#####################################

# Here we compute the overall trends and number of misses outside of the bootstrap


### Compute visit counts ------------------------------------------------------
sim_tm_all <- all_dx_visits %>% 
  mutate(period = -days_since_index) %>%
  distinct(enrolid,period,days_since_index) %>%
  inner_join(sim_obs,by = c("enrolid", "days_since_index")) 

all_vis_count <- sim_tm_all %>% 
  count(.,period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

sim_tm_ssd <- all_dx_visits %>% 
  mutate(period = -days_since_index) %>%
  inner_join(ssd_codes,by = c("dx", "dx_ver")) %>% 
  distinct(enrolid,period,days_since_index) %>%
  inner_join(sim_obs,by = c("enrolid", "days_since_index")) 

ssd_vis_count <- sim_tm_ssd %>% 
  count(.,period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))


#### fit trends ----------------------------------------------------------------

all_vis_count <- return_fits(all_vis_count,
                             model = delay_params$final_model,
                             cp = delay_params$cp,
                             periodicity = delay_params$periodicity)

ssd_vis_count <- return_fits(ssd_vis_count,
                             model = delay_params$final_model,
                             cp = delay_params$cp,
                             periodicity = delay_params$periodicity)


all_vis_count %>% summarise(n_miss = sum(num_miss,na.rm = T))
ssd_vis_count %>% summarise(n_miss = sum(num_miss,na.rm = T))


################################
#### Prepare Bootstrap Data ####
################################

#### Setup outer bootstrap data ------------------------------------------------

# sample overall set of patients
tmp <- tibble(boot_trial=1:delay_params$boot_trials) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(enrolid = sample(patient_ids$enrolid, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

#### compute counts for all visits ---------------------------------------------
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "enrolid") %>%
                        mutate(period = -days_since_index) %>%
                        distinct(enrolid,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("enrolid", "days_since_index")) %>% 
                        rename(patient_id=enrolid)))

tmp <- tmp %>% 
  mutate(all_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7))))

#### compute counts for ssd visits ---------------------------------------------
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "enrolid") %>%
                        mutate(period = -days_since_index) %>%
                        inner_join(ssd_codes,by = c("dx", "dx_ver")) %>%
                        distinct(enrolid,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("enrolid", "days_since_index")) %>% 
                        rename(patient_id=enrolid)))

tmp <- tmp %>% 
  mutate(ssd_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7))))

# remove the sim_tm data
tmp <- tmp %>% 
  select(boot_trial,boot_sample,all_vis_count,ssd_vis_count)


### fit models -----------------------------------------------------------------

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


boot_data <- tmp
rm(tmp)

# save bootstrap data
save(boot_data,file = paste0(sim_out_path,"boot_data.RData"))

#### Compute number missed -----------------------------------------------------
boot_data %>% 
  summarise(mean_n_miss_all = mean(n_miss_all),
            median_n_miss_all = median(n_miss_all),
            mean_n_miss_ssd = mean(n_miss_ssd),
            median_n_miss_ssd = median(n_miss_ssd))



################################
#### Run bootstrap analysis ####
################################


cluster <- parallel::makeCluster(33)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterCall(cluster, function() library(delaySim))
parallel::clusterExport(cluster,c("run_boot_trials","boot_data","sim_tm_all","sim_tm_ssd","delay_params","n_patients"),
                        envir=environment())


sim_res_all <- parallel::parLapply(cl = cluster,
                               1:delay_params$boot_trials,
                               function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                           tm_data = sim_tm_all,
                                                           miss_bins = select(boot_data$all_vis_count[[x]],period,num_miss),
                                                           delay_params = delay_params,
                                                           n_patients = n_patients)})

sim_res_all <- map2(sim_res_all,1:delay_params$boot_trials,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

sim_res_ssd <- parallel::parLapply(cl = cluster,
                                   1:delay_params$boot_trials,
                                   function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                               tm_data = sim_tm_ssd,
                                                               miss_bins = select(boot_data$ssd_vis_count[[x]],period,num_miss),
                                                               delay_params = delay_params,
                                                               n_patients = n_patients)})

sim_res_ssd <- map2(sim_res_ssd,1:delay_params$boot_trials,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

parallel::stopCluster(cluster)

rm(cluster)

save(sim_res_all,file = paste0(sim_out_path,"sim_res_all.RData"))
save(sim_res_ssd,file = paste0(sim_out_path,"sim_res_ssd.RData"))

# now extract the specific sim obs values for future steps
tmp1 <- sim_res_ssd %>% unnest(res) %>% distinct(obs)
tmp2 <- sim_res_all %>% unnest(res) %>% distinct(obs)

sim_obs_reduced <- sim_obs %>% 
  inner_join(distinct(bind_rows(tmp1,tmp2))) %>% 
  mutate(period = -days_since_index) %>% 
  select(obs,enrolid,period)

save(sim_obs_reduced, file = paste0(sim_out_path,"sim_obs_reduced.RData"))

# cleanup
rm(list = ls()[!(ls() %in% c("delay_params","sim_out_path","sim_in_path","cond_name"))])
gc()

###################################
#### Analyze Bootstrap Results ####
###################################

# Load in sumulation results
load(paste0(sim_out_path,"sim_res_ssd.RData"))
load(paste0(sim_out_path,"sim_res_all.RData"))
load(paste0(sim_out_path,"sim_obs_reduced.RData"))

sim_res_ssd$res[[1]] %>% 
  inner_join(sim_obs_reduced) %>% 
  group_by(enrolid,boot_id) %>% 
  summarise(duration=max(period),
            n_miss=n())

options(dplyr.summarise.inform = FALSE)

compute_boot_stats <- function(sim_data,sim_obs_data,delay_params,n_patients){
  
  tmp1 <- sim_data %>% 
    inner_join(sim_obs_data, by = "obs") %>% 
    group_by(enrolid,boot_id) %>%
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

compute_boot_stats(sim_res_ssd$res[[1]], sim_obs_reduced, delay_params = delay_params, n_patients = n_patients)



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


########################
#### Export Results ####
########################

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

aggregate_sim_stats(sim_stats_ssd)
aggregate_sim_stats(sim_stats_all)


sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(mean(.),2))


sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2))

sim_stats_ssd %>% 
  select(sim_trial,boot_trial,miss_tab) %>% 
  unnest(miss_tab) %>% 
  group_by(miss_bin) %>% 
  summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2))
