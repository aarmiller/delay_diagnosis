


rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")

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


########################################
###### Apply Filtering Criteria ########
########################################

# Apply any criteria to restrict to final study population here


######################
#### Overall Fits ####
######################


#### Load Data  ----------------------------------------------------------------

load(paste0(sim_in_path,"all_dx_visits.RData"))
load(paste0(sim_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver)


patient_ids <- tm %>% 
  distinct(enrolid)

n_patients <- nrow(patient_ids)

#### Compute visit counts ------------------------------------------------------
sim_tm_all <- all_dx_visits %>% 
  mutate(period = -days_since_index) %>%
  distinct(enrolid,period,days_since_index) %>%
  inner_join(sim_obs,by = c("enrolid", "days_since_index")) %>% 
  rename(patient_id=enrolid)

all_vis_count <- sim_tm_all %>% 
  count(.,period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

sim_tm_ssd <- all_dx_visits %>% 
  mutate(period = -days_since_index) %>%
  inner_join(ssd_codes,by = c("dx", "dx_ver")) %>% 
  distinct(enrolid,period,days_since_index) %>%
  inner_join(sim_obs,by = c("enrolid", "days_since_index")) %>% 
  rename(patient_id=enrolid)

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

tmp <- tibble(boot_trial=1:delay_params$boot_trials) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(enrolid = sample(patient_ids$enrolid, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# compute counts for all visits
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

# compute counts for ssd visits
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

tmp <- tmp %>% 
  select(boot_trial,boot_sample,all_vis_count,ssd_vis_count)

tmp$all_vis_count[[1]]
tmp$ssd_vis_count[[1]]


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
                          ~summarise(.,n_miss_all = sum(num_miss,na.rm = T)))) %>% 
  unnest(n_miss_all) %>% 
  mutate(n_miss_ssd = map(ssd_vis_count,
                          ~summarise(.,n_miss_ssd = sum(num_miss,na.rm = T)))) %>% 
  unnest(n_miss_ssd)


tmp %>% 
  summarise(n_miss_all = mean(n_miss_all),
            n_miss_ssd = mean(n_miss_ssd))



################################
#### Run bootstrap analysis ####
################################

runs_sim_draw_visits()

sim_tm_all

select(tmp$all_vis_count[[1]],period,num_miss)


tmp_sim_data <- list(time_map = mutate(sim_tm_all,miss_ind=1),
                     miss_bins_visits = select(tmp$all_vis_count[[1]],period,num_miss),
                     change_point = delay_params$cp,
                     total_patients = n_patients,
                     miss_bins = delay_params$miss_bins,
                     duration_bins = delay_params$duration_bins)

run_sim_draw_simple(tmp_sim_data)

mutate(sim_tm_all,miss_ind=1) %>% 
  dplyr::filter(miss_ind == 1) %>% 
  dplyr::inner_join(select(tmp$all_vis_count[[1]],period,num_miss), by = "period") %>% 
  dplyr::mutate(rand = runif(n = dplyr::n())) %>% 
  dplyr::arrange(period, rand) %>% 
  dplyr::group_by(period) %>% 
  dplyr::filter(dplyr::row_number() <= num_miss) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-num_miss, -rand)

select(tmp$all_vis_count[[1]],period,num_miss) %>% summarise(sum(num_miss))
tmp

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






return_fits <- function(data, model, cp, periodicity = TRUE){
  
  if (periodicity==TRUE){
    
    if (model == "linear"){
      
      fit <- lm(n ~ period + dow, filter(data, period>cp))
      
    } else if (model == "exponential") {
      
      fit <- lm(n ~ log(period) + dow, filter(data, period>cp))
      
    } else if (model == "quadratic") {
      
      fit <- lm(n ~ poly(period,degree = 2) + dow, filter(data, period>cp))
      
    } else if (model == "cubic") {
      
      fit <- lm(n ~ poly(period,degree = 3) + dow, filter(data, period>cp))
      
    }
    
  } else {
    
    if (model == "linear"){
      
      fit <- lm(n ~ period, filter(data, period>cp))
      
    } else if (model == "exponential") {
      
      fit <- lm(n ~ log(period), filter(data, period>cp))
      
    } else if (model == "quadratic") {
      
      fit <- lm(n ~ poly(period,degree = 2), filter(data, period>cp))
      
    } else if (model == "cubic") {
      
      fit <- lm(n ~ poly(period,degree = 3), filter(data, period>cp))
      
    }
    
  }
  
  data %>%
    mutate(pred = predict(fit,newdata = .)) %>% 
    mutate(num_miss = n-pred) %>%
    mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
    mutate(num_miss = ifelse(period>cp,NA,num_miss))
  
}
  

#### fit models ----------------------------------------------------------------
if (delay_params$periodicity){
  # if periodicity = TRUE add day of week term
  fit1 <- lm(n ~ period + dow, filter(all_vis_count, period>delay_params$cp))
  fit2 <- lm(n ~ log(period) + dow, filter(all_vis_count, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2) + dow, filter(all_vis_count, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3) + dow, filter(all_vis_count, period>delay_params$cp))
} else {
  fit1 <- lm(n ~ period, filter(all_vis_count, period>delay_params$cp))
  fit2 <- lm(n ~ log(period), filter(all_vis_count, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2), filter(all_vis_count, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3), filter(all_vis_count, period>delay_params$cp))
}