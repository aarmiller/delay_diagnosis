
rm(list =ls())
library(tidyverse)
library(bit64)


cond_name <- "sepsis_pre_covid"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)

# final models that were run
models <- tibble(cp = delay_params$cp) %>% 
  mutate(model = map(cp,~delay_params$final_model)) %>% 
  unnest(model)

sim_out_path <- paste0(delay_params$out_path,"sim_results/")

# load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

# extract patient ids and number of patients
patient_ids <- index_cases %>% 
  distinct(patient_id)

n_patients <- nrow(patient_ids)


# Report data path

report_data_path <- paste0(sim_out_path,"report_data/")

dir.create(report_data_path)


#######################
#### Baseline Data ####
#######################


load(paste0(delay_params$base_path,"delay_results/demo_data.RData"))

age_cats <- c(-1,17,34,44,54,64,130)

tmp1 <- demo1 %>% 
  inner_join(index_cases) %>% 
  count(sex) %>% 
  mutate(pct = round(100*(n/n_patients),1)) %>% 
  mutate(sex = ifelse(sex==1, "Male","Female")) %>% 
  rename(variable = sex)

tmp2 <- demo1 %>% 
  inner_join(index_cases) %>% 
  mutate(age = year(as_date(index_date))-dobyr) %>% 
  mutate(age_cat = as.character(cut(age,age_cats))) %>% 
  count(age_cat) %>% 
  mutate(pct = round(100*(n/n_patients),1)) %>% 
  inner_join(tribble(~age_cat,~variable,
                     '(-1,17]', 'Age <18',
                     '(17,34]', 'Age 18-34',
                     '(34,44]', 'Age 35-44',
                     '(44,54]', 'Age 45-54',
                     '(54,64]', 'Age 55-64',
                     '(64,130]', 'Age >= 65')) %>% 
  select(variable,n,pct)

tmp3 <- demo2 %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  inner_join(distinct(index_cases,patient_id)) %>% 
  count(source) %>% 
  mutate(pct = round(100*(n/n_patients),1)) %>% 
  inner_join(tribble(~source,~variable,
                     "ccae", "Commercial Insurance",
                     "mdcr", "Medicare",
                     "medicaid", "Medicaid")) %>% 
  select(variable,n,pct)

tmp4 <- index_cases %>% 
  mutate(dow=weekdays(as_date(index_date))) %>% 
  mutate(dow = fct_relevel(dow,"Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")) %>% 
  count(dow) %>% 
  mutate(pct = round(100*(n/n_patients),1)) %>% 
  mutate(variable = as.character(dow)) %>% 
  select(variable,n,pct)

tmp5 <- index_cases %>% 
  mutate(year = year(as_date(index_date))) %>% 
  count(year) %>% 
  mutate(pct = round(100*(n/n_patients),1)) %>% 
  mutate(variable=paste0("Year ",year)) %>% 
  select(variable,n,pct)


pat_char <- bind_rows(tmp1,tmp2,tmp3)

index_char <- bind_rows(tmp4,tmp5)

baseline_data <- list(pat_char = pat_char,
                      index_char = index_char)

baseline_data$index_char

save(baseline_data, file = paste0(report_data_path,"baseline_data.RData"))
  

# load time map
load(paste0(delay_params$base_path,"delay_results/delay_tm.RData"))
tm <- inner_join(tm, patient_ids, by = join_by(patient_id))

tm %>% 
  filter(days_since_index==0)


################
#### Remove ####
################

# run markdown report
rmarkdown::render(input = "github/delay_diagnosis/projects/sepsis_revised10/pre_covid/final_delay_report_sepsis_pre_covid.Rmd",
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/"))


######################
#### Trends Info #####
######################
# 
load(paste0(sim_out_path,"trends/fit_trends.RData"))


trends_ssd <- ssd_vis_count_fitted %>%
  rename(tot_miss=num_miss) %>%
  unnest(counts) %>%
  mutate(model_label = paste0(model," (CP = ",cp," days)"))

trends_all <- all_vis_count_fitted %>%
  rename(tot_miss=num_miss) %>%
  unnest(counts) %>%
  mutate(model_label = paste0(model," (CP = ",cp," days)"))

model_comparisons_table <- bind_rows(trends_ssd %>% 
            group_by(model,cp) %>% 
            filter(is.na(num_miss)) %>% 
            summarise(`RMSE (before CP)`=sqrt(mean((n-pred)^2)),
                      `Total Misses` = mean(tot_miss)) %>% 
            mutate(group = "SSD Visits"),
          trends_all %>% 
            group_by(model,cp) %>% 
            filter(is.na(num_miss)) %>% 
            summarise(`RMSE (before CP)`=sqrt(mean((n-pred)^2)),
                      `Total Misses` = mean(tot_miss)) %>% 
            mutate(group = "All Visits")) %>% 
  select(group,everything()) %>% ungroup()


trends_out <- list(model_comparisons_table = model_comparisons_table)

save(trends_out, file = paste0(report_data_path,"trends_output.RData"))

##########################################
#### Aggregated Simulation Statistics ####
##########################################



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

##########################
#### Loop over models ####
##########################

load(paste0(delay_params$base_path,"delay_results/all_dx_visits.RData"))

folder_list <- c("exponential_cp7","exponential_cp14")

i <- "exponential_cp14"

for (i in folder_list){
  tmp_in_path <- paste0(sim_out_path,i,"/")
  
  load(paste0(tmp_in_path,"sim_res_ssd.RData"))
  load(paste0(tmp_in_path,"sim_res_all.RData"))
  
  ## filter to sim obs used for quicker merge ------------------------------------
  
  tmp1 <- sim_res_ssd %>% 
    mutate(obs = map(res,~distinct(.,obs))) %>% 
    select(obs) %>% 
    unnest(obs) %>% 
    distinct(obs)
  
  tmp2 <- sim_res_all %>% 
    mutate(obs = map(res,~distinct(.,obs))) %>% 
    select(obs) %>% 
    unnest(obs) %>% 
    distinct(obs)
  
  tmp <- bind_rows(tmp1,tmp2) %>% distinct(obs)
  
  sim_obs_reduced <- inner_join(tmp,sim_obs) %>% 
    mutate(period = -days_since_index) %>% 
    select(-days_since_index)
  
  rm(tmp1,tmp2,tmp)
  gc()
  
  ### Aggregate results --------------------------------------------------------
  # setup cluster
  cluster <- parallel::makeCluster(30)
  parallel::clusterCall(cluster, function() library(tidyverse))
  parallel::clusterExport(cluster,c("compute_boot_stats","sim_obs_reduced","delay_params","n_patients"),
                          envir=environment())
  
  tmp1 <- parallel::parLapply(cl = cluster,
                              sim_res_ssd$res,
                                     function(x){compute_boot_stats(x,
                                                                    sim_obs_reduced,
                                                                    delay_params = delay_params,
                                                                    n_patients = n_patients)})
  
  sim_stats_ssd <- tmp1 %>% enframe() %>%
    mutate(stats = map(value,~.$stats)) %>%
    mutate(miss_tab = map(value,~.$miss_tab)) %>%
    mutate(dur_tab = map(value,~.$dur_tab)) %>% 
    select(-value)
  
  rm(tmp1,sim_res_ssd)
  gc()
  
  # aggregate All visit results
  tmp2 <- parallel::parLapply(cl = cluster,
                              sim_res_all$res,
                              function(x){compute_boot_stats(x,
                                                             sim_obs_reduced,
                                                             delay_params = delay_params,
                                                             n_patients = n_patients)})
  
  sim_stats_all <- enframe(tmp2) %>%
    mutate(stats = map(value,~.$stats)) %>%
    mutate(miss_tab = map(value,~.$miss_tab)) %>%
    mutate(dur_tab = map(value,~.$dur_tab)) %>% 
    select(-value)
  
  save(sim_stats_ssd,sim_stats_all,file = paste0(tmp_in_path,"sim_stats.RData"))
  
  rm(tmp2,sim_stats_all,sim_stats_ssd,sim_res_all)
  gc()
  
}


###################################
#### Export Simulation Results ####
###################################



load(paste0(delay_params$out_path,"index_cases.RData"))

n_patients <- nrow(index_cases)

### Generate aggregate statistics ----------------------------------------------

aggregate_sim_stats <- function(sim_stats_data, n_patients){
  
  ## Aggregate summary stats ---------------------------------------------------
  
  tmp1 <- sim_stats_data %>%
    select(stats) %>%
    unnest(stats) %>%
    summarise_at(vars(n_pat:median_n_miss),~round(mean(.),2))

  tmp2 <- sim_stats_data %>%
    select(stats) %>%
    unnest(stats) %>%
    summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.025)),2))

  tmp3 <- sim_stats_data %>%
    select(stats) %>%
    unnest(stats) %>%
    summarise_at(vars(n_pat:median_n_miss),~round(quantile(.,probs = c(0.975)),2))

  res1 <- inner_join(gather(tmp1, key=measure, value = mean),
                     gather(tmp2, key=measure, value = low),
                     by = join_by(measure)) %>%
    inner_join(gather(tmp3, key=measure, value = high),
               by = join_by(measure)) 
  
  
  res2 <- res1 %>% 
    filter(measure=="n_pat") %>% 
    mutate_at(vars(mean:high),~round(100*./n_patients,2)) %>% 
    mutate(measure_out=paste0(mean," (",low,"-",high,")")) %>% 
    mutate(measure="pct_pat")
  
  out1 <- bind_rows(res2,res1) %>% 
    mutate(measure_out=paste0(mean," (",low,"-",high,")"))
  
  ## Aggregate miss bins -------------------------------------------------------
  
  tmp1 <- sim_stats_data %>% 
    select(miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(mean(.),2))
  
  tmp2 <- sim_stats_data %>% 
    select(miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2))
  
  tmp3 <- sim_stats_data %>% 
    select(miss_tab) %>% 
    unnest(miss_tab) %>% 
    group_by(miss_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2))
  
  out2 <- tmp1 %>% 
    gather(key = key, value = mean,-miss_bin) %>% 
    inner_join(gather(tmp2, key = key, value = low,-miss_bin), by = join_by(miss_bin, key)) %>% 
    inner_join(gather(tmp3, key = key, value = high,-miss_bin), by = join_by(miss_bin, key)) %>% 
    mutate(measure_out=paste0(mean," (",low,"-",high,")")) %>% 
    select(miss_bin,key,measure_out) %>% 
    spread(key = key,value = measure_out)
  
  ## Aggregate duration bins ---------------------------------------------------
  
  tmp1 <- sim_stats_data %>% 
    select(dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(mean(.),2))
  
  tmp2 <- sim_stats_data %>% 
    select(dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2))
  
  tmp3 <- sim_stats_data %>% 
    select(dur_tab) %>% 
    unnest(dur_tab) %>% 
    group_by(duration_bin) %>% 
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2))
  
  out3 <- tmp1 %>% 
    gather(key = key, value = mean,-duration_bin) %>% 
    inner_join(gather(tmp2, key = key, value = low,-duration_bin), by = join_by(duration_bin, key)) %>% 
    inner_join(gather(tmp3, key = key, value = high,-duration_bin), by = join_by(duration_bin, key)) %>% 
    mutate(measure_out=paste0(mean," (",low,"-",high,")")) %>% 
    select(duration_bin,key,measure_out) %>% 
    spread(key = key,value = measure_out)
  

  return(list(main_stats = out1,
              miss_bins = out2,
              duration_bins = out3))
}

load(paste0(sim_out_path,"exponential_cp14/sim_stats.RData"))
agg_stats_ssd_exp14 <- aggregate_sim_stats(sim_stats_ssd,n_patients)
agg_stats_all_exp14 <- aggregate_sim_stats(sim_stats_all,n_patients)

load(paste0(sim_out_path,"exponential_cp7/sim_stats.RData"))
agg_stats_ssd_exp7 <- aggregate_sim_stats(sim_stats_ssd,n_patients)
agg_stats_all_exp7 <- aggregate_sim_stats(sim_stats_all,n_patients)


agg_stats <- list(exponential_7 = list(agg_stats_ssd = agg_stats_ssd_exp7,
                                        agg_stats_all = agg_stats_all_exp7),
                  exponential_14 = list(agg_stats_ssd = agg_stats_ssd_exp14,
                                         agg_stats_all = agg_stats_all_exp14))

save(agg_stats, file = paste0(report_data_path, "agg_stats.RData"))


