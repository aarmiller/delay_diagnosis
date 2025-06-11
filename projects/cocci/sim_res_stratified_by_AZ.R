library(tidyverse)
library(lubridate)
library(smallDB)
library(furrr)
library(bit64)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "cocci"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 91 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Connect to db
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

### Load index cases
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),"index_cases.RData"))

num_cores <- 50
cp_selected <- delay_params$cp - 1 # minus 1 as the risk factors are for within delay window

#update demo1 and demo2
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
demo2 <- demo2 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")

rural_ids <- rural_visits %>% inner_join(index_cases %>% distinct(patient_id)) %>% distinct(patient_id)

# update all_dx_visits
load(paste0(delay_base_path,"all_dx_visits.RData"))
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

# update sim_obs
sim_obs <- sim_obs %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<0)

# update time map
load(paste0(delay_base_path,"delay_tm.RData"))

problem_patient_ids <- tm %>% filter(days_since_index == 0) %>% distinct(patient_id) %>%
  anti_join(tm %>% distinct(patient_id),.)

tm <- tm %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

#Load in the setting types
load(paste0(delay_base_path,"visit_info.RData"))

# update tm_stdplac
tm_stdprov <- tm_stdprov %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

load(paste0(sim_out_path,"/sim_res_ssd.RData"))
sim_res_ssd <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  select(trial, boot_trial, res) %>% 
  unnest(res) %>% 
  select(-boot_id)

load(paste0(sim_out_path,"sim_obs_reduced.RData"))

sim_res_sim_obs <- sim_res_ssd %>% 
  inner_join(sim_obs_reduced, by = "obs") 

load(paste0(sim_out_path,"boot_data.RData"))
boot_data <- boot_data %>% select(boot_trial, boot_sample)

n_trials <- nrow(distinct(sim_res_ssd,trial))

# Demo data
reg_demo <- demo1 %>% 
  mutate(female=(sex==2),
         age = index_year-dobyr,
         stdrace = as.numeric(stdrace)) %>% 
  left_join(tibble(stdrace = c(0,1,2,4,9),
                   race = c("Missing/Unknown","White","Black","Hispanic","Other")),
            by = "stdrace") %>% 
  mutate(race = fct_relevel(race,"White"))

reg_demo <- reg_demo %>% 
  left_join(demo2 %>% 
              filter(index_date<=dtend & index_date>=dtstart) %>% 
              mutate(msa_new = msa %in% c("0","")) %>% 
              mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
              mutate(source = as.factor(source)) %>% 
              distinct(patient_id,source,msa=msa_new),
            by = "patient_id")

age_cats <- c(-1,17,34,44,54,64,130)

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,age,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))

tm <- tm %>% rename(admdate = svcdate)

## Prepare prior geography  ----------------------------------------------------
# source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
# load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024

source(paste0(stringr::str_replace(delay_params$out_path,
                                   paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
              "get_enroll_detail_fun.R"))
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
            "egeoloc_labels.RData"))

enroll_collapsed_temp <- gather_collapse_enrollment(enrolid_list = reg_demo %>% distinct(patient_id) %>% .$patient_id,
                                                    vars = "egeoloc",
                                                    db_path =  paste0(delay_params$small_db_path,"cocci.db"),
                                                    num_cores=10,
                                                    collect_tab = collect_table(year = 1:22))

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322071/
enroll_collapsed_temp2 <- enroll_collapsed_temp %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  distinct() 
# only 26,378 out of 26,905 have location information 526 are medicaid 1 missing

AZ_location_ind <- enroll_collapsed_temp2 %>%  
  mutate(AZ = ifelse(state_abb %in% c("AZ"), 1L, ifelse(is.na(state_abb), NA, 0))) %>% 
  select(patient_id, location:state_abb, AZ) %>% 
  distinct()
# location_ind %>% filter(top2_high_inc_state_baddley == 1) %>% distinct(state_abb) #check coding
# 26,048 of the 26,378 with location info have non missing state, 330 missing state obs

reg_demo <- left_join(reg_demo, AZ_location_ind)

reg_demo %>% count(AZ, location) %>% 
  print(n = Inf) # notice that when location = puerto rico, AZ is NA to lets fix that

reg_demo <- reg_demo %>% mutate(AZ = ifelse(is.na(AZ) & location == "puerto rico", 0L, AZ))

reg_demo %>% count(AZ, location) %>% 
  print(n = Inf)  # now when location = puerto rico, AZ is 0

reg_demo %>% count(AZ, location) %>% 
  write_csv(paste0(sim_out_path,"AZ_ind_counts.csv")) # Here the 527 NA NA are the 526 medicaid and 1 no location obs

reg_demo %>% count(AZ, source) # 833 total patients without state info of which 526 are medicaid,
# 212 are ccae, 95 are mdcr

reg_demo %>% filter(AZ == 1) %>% 
  count(AZ, location, state_name,state_abb) #check

reg_demo %>% filter(AZ ==0 | is.na(AZ)) %>% 
  count(AZ, location, state_name,state_abb) %>%  #check
  print(n = Inf)

AZ_ind_data <- reg_demo
save(AZ_ind_data, file = paste0(sim_out_path,"AZ_ind_data.RData"))

# Subset to project specific patient ids  ######################################
index_cases_main <- index_cases

for (x_val in unique(reg_demo$AZ)){

  rm(list = ls()[!ls() %in% c("index_cases_main", "x_val", "reg_demo",
                              "delay_params", "sim_out_path")])
  
  load(paste0(sim_out_path,"sim_res_ssd.RData"))
  load(paste0(sim_out_path,"sim_obs_reduced.RData"))
  load(paste0(delay_params$base_path,"/delay_results/all_dx_visits.RData"))
  load(paste0(delay_params$base_path,"/delay_results/delay_tm.RData"))
  load(paste0(stringr::str_replace(delay_params$out_path,
                                   paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
              "index_cases.RData"))
  load(paste0(sim_out_path,"boot_data.RData"))
  
  # Subset to project specific patient ids
  if(is.na(x_val)) {
    index_cases <- index_cases_main %>% inner_join(reg_demo %>% filter(is.na(AZ)) %>% distinct(patient_id) )
  } else {
    index_cases <- index_cases_main %>% inner_join(reg_demo %>% filter(AZ == x_val) %>% distinct(patient_id) )
  }
 
  index_dx_dates <- index_cases
  patient_ids <- index_cases %>% distinct(patient_id)
  n_patients <- nrow(patient_ids)
  
  # update all_dx_visits
  all_dx_visits <- all_dx_visits %>%
    inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
    mutate(days_since_index = days_since_index - shift) %>% 
    select(-shift) %>% 
    filter(days_since_index<=0) 
  
  # update time map
  tm <- tm %>% 
    inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
    mutate(days_since_index = days_since_index - shift) %>% 
    select(-shift) %>%
    filter(days_since_index<=0)
  
  sim_obs_reduced <- sim_obs_reduced %>% 
    inner_join(index_cases %>% select(patient_id), by = "patient_id") 
  
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
  plan(multisession, workers = 40)
  
  tmp1 <- sim_res_ssd %>% 
    # slice(1:30) %>% 
    mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
                                                           sim_obs_reduced,
                                                           delay_params = delay_params, 
                                                           n_patients = n_patients)))
  plan(sequential)
  
  tmp1 <- tmp1 %>% 
    select(sim_trial,boot_trial,sim_stats) %>% 
    mutate(stats = map(sim_stats,~.$stats)) %>% 
    mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>% 
    mutate(dur_tab = map(sim_stats,~.$dur_tab))
  
  sim_stats_ssd <- tmp1 %>% 
    select(-sim_stats)
  
  save(sim_stats_ssd,file = paste0(sim_out_path,"sim_stats_AZ_", x_val , ".RData"))
  
  
  ###################################
  #### Export Simulation Results ####
  ###################################
  
  # Note: We may want to move parts of this section to a separate markdown report script
  
  ### Generate aggregate statistics ----------------------------------------------
  
  aggregate_sim_stats <- function(sim_stats_data){
    tmp1 <- sim_stats_data %>% 
      select(sim_trial,boot_trial,stats) %>% 
      unnest(stats) %>% 
      summarise_at(vars(n_pat:median_n_miss),~round(mean(., na.rm = T),2))
    
    tmp2 <- sim_stats_data %>% 
      select(sim_trial,boot_trial,stats) %>% 
      unnest(stats) %>% 
      summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.025), na.rm = T),2))
    
    tmp3 <- sim_stats_data %>% 
      select(sim_trial,boot_trial,stats) %>% 
      unnest(stats) %>% 
      summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.975), na.rm = T),2))
    
    
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
  
  ### Generate bins for number of misses per patient -----------------------------
  
  # SSD Visits
  tmp1 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(mean(., na.rm = T),2)) %>% 
    gather(key = key, value = mean, -miss_bin) %>% 
    mutate(mean = round(mean,2))
  
  tmp2 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(.,, na.rm = T, probs = c(0.025)),2)) %>% 
    gather(key = key, value = low, -miss_bin) %>% 
    mutate(low = round(low,2))
  
  tmp3 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(., na.rm = T, probs = c(0.975)),2)) %>% 
    gather(key = key, value = high, -miss_bin) %>% 
    mutate(high = round(high,2))
  
  miss_bins_ssd <- inner_join(tmp1,tmp2) %>% 
    inner_join(tmp3) %>% 
    mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
    select(miss_bin,key,out) %>% 
    spread(key = key, value = out)
  
  
  ### Generate bins for duration of misses ---------------------------------------
  
  # SSD Visits
  tmp1 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(mean(., na.rm = T),2)) %>% 
    gather(key = key, value = mean, -duration_bin) %>% 
    mutate(mean = round(mean,2))
  
  tmp2 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(., na.rm = T, probs = c(0.025)),2)) %>% 
    gather(key = key, value = low, -duration_bin) %>% 
    mutate(low = round(low,2))
  
  tmp3 <- sim_stats_ssd %>% 
    select(sim_trial,boot_trial,dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(., na.rm = T, probs = c(0.975)),2)) %>% 
    gather(key = key, value = high, -duration_bin) %>% 
    mutate(high = round(high,2))
  
  dur_bins_ssd <- inner_join(tmp1,tmp2) %>% 
    inner_join(tmp3) %>% 
    mutate(out = paste0(mean," (",low,"-",high,")")) %>% 
    select(duration_bin,key,out) %>% 
    spread(key = key, value = out)
  
  
  save(agg_stats_ssd,
       miss_bins_ssd,
       dur_bins_ssd, 
       file = paste0(sim_out_path,"aggregated_sim_results_AZ_", x_val , ".RData"))
}

