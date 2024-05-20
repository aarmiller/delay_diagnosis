
library(tidyverse)
library(smallDB)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "histo"

con <- DBI::dbConnect(RSQLite::SQLite(), paste0("/Shared/AML/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

delay_base_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")

load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]
rm(delay_any_params)

age_cats <- c(-1,17,34,44,54,64,130)

##################
### Load Data ####
##################

load(paste0(delay_base_path,"delay_tm.RData"))
load(paste0(delay_base_path,"demo_data.RData"))
load(paste0(delay_base_path,"ssd_visit/sim_res.RData"))

abx_codes <- unlist(codeBuildr::load_rx_codes("all_abx"),use.names = F)

index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=delay_params$upper_bound)


abx_vis <- con %>% 
  tbl("all_rx_visits") %>% 
  collect() %>% 
  filter(ndcnum %in% abx_codes)



##################################
#### Assemble Data for Models ####
##################################

# Compute duration
sim_dur <- sim_res_sim_obs %>% 
  inner_join(sim_res) %>% 
  group_by(patient_id,trial) %>% 
  summarise(duration = -min(days_since_index)) %>% 
  ungroup()

reg_data_base <- demo1 %>% 
  inner_join(index_dx_dates) %>% 
  mutate(age = index_year-dobyr) %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  select(patient_id,sex,age_cat,index_year,index_date) %>% 
  mutate(sex = ifelse(sex==2,"female","male"))



###############################
#### Any Antibiotic Models ####
###############################

abx_inds <- reg_data_base %>% 
  inner_join(distinct(abx_vis,patient_id,date)) %>% 
  filter(date<index_date,date>=(index_date-delay_params$cp)) %>% 
  mutate(time_to_index = index_date-date) %>% 
  distinct(patient_id,time_to_index) %>% 
  mutate(abx_vis=1L)


# restrict to only antibiotics during individual delay window
run_dur_model1 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num)
  
  tmp_inds <- abx_inds %>% 
    inner_join(tmp,by = join_by(patient_id)) %>% 
    filter(time_to_index<=duration) %>% 
    distinct(patient_id,abx_vis)
  
  tmp_reg_data <- tmp %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_inds,by = join_by(patient_id)) %>% 
    mutate(abx_vis=replace_na(abx_vis,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_vis, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
  
}
# run_dur_model1(4)

# Don't restrict to individual delay window, use entire window
run_dur_model1_2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num)
  
  tmp_inds <- abx_inds %>% 
    inner_join(tmp,by = join_by(patient_id)) %>% 
    # filter(time_to_index<=duration) %>% 
    distinct(patient_id,abx_vis)
  
  tmp_reg_data <- tmp %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_inds,by = join_by(patient_id)) %>% 
    mutate(abx_vis=replace_na(abx_vis,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_vis, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}
# run_dur_model1_2(3)

## Run Model 1 -----------------------------------------------------------------
mod_res1 <- distinct(sim_res,trial)

mod_res1 <- mod_res1 %>% 
  mutate(res = map(trial,run_dur_model1))

## Run Model 2 -----------------------------------------------------------------
mod_res1_2 <- distinct(sim_res,trial)

mod_res1_2 <- mod_res1 %>% 
  mutate(res = map(trial,run_dur_model1_2))

## Assemble Results ------------------------------------------------------------

mod_res1 <- mod_res1 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mdn_est = median(value),
            mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

mod_res1_2 <- mod_res1_2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mdn_est = median(value),
            mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

#####################################
#### Days of antibiotic duration ####
#####################################

abx_days <- abx_vis %>% 
  filter(daysupp>0) %>% 
  inner_join(select(reg_data_base,patient_id,index_date)) %>% 
  filter(date >= index_date-delay_params$cp,
         date < index_date) %>% 
  select(patient_id:daysupp)

### Functions ------------------------------------------------------------------

# using max duration of antibiotic recieved during individual diagnostic opportunity window
run_dur_model2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num) %>% 
    inner_join(select(reg_data_base,patient_id,index_date),by = join_by(patient_id))
  
  tmp_abx_days <- abx_days %>% 
    inner_join(tmp, by = join_by(patient_id)) %>% 
    filter(date<index_date,
           date>=(index_date-duration)) %>% 
    group_by(patient_id) %>% 
    summarise(daysupp = max(daysupp)) %>% 
    rename(abx_days=daysupp)
  
  tmp_reg_data <- tmp %>% 
    select(patient_id,duration) %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_abx_days,by = join_by(patient_id)) %>% 
    mutate(abx_days=replace_na(abx_days,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_days, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}
# run_dur_model2(4)

# using max duration of antibiotic recieved during individual diagnostic opportunity window
run_dur_model2_2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num) %>% 
    inner_join(select(reg_data_base,patient_id,index_date),by = join_by(patient_id))
  
  tmp_abx_days <- abx_days %>% 
    inner_join(tmp, by = join_by(patient_id)) %>% 
    group_by(patient_id) %>% 
    summarise(daysupp = max(daysupp)) %>% 
    rename(abx_days=daysupp)
  
  tmp_reg_data <- tmp %>% 
    select(patient_id,duration) %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_abx_days,by = join_by(patient_id)) %>% 
    mutate(abx_days=replace_na(abx_days,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_days, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}
# run_dur_model2_2(4)

### Run Model 1 ----------------------------------------------------------------

mod_res2 <- distinct(sim_res,trial)

mod_res2 <- mod_res2 %>% 
  mutate(res = map(trial,run_dur_model2))

### Run Model 2 ----------------------------------------------------------------

mod_res2_2 <- distinct(sim_res,trial)

mod_res2_2 <- mod_res2_2 %>% 
  mutate(res = map(trial,run_dur_model2_2))

## Assemble Results ------------------------------------------------------------

mod_res2 <- mod_res2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mdn_est = median(value),
            mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

mod_res2_2 <- mod_res2_2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mdn_est = median(value),
            mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))



###########################################
#### Categories of antibiotic duration ####
###########################################


mod_res %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))


## binning duration ------------------------------------------------------------


run_dur_model3 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num) %>% 
    inner_join(select(reg_data_base,patient_id,index_date),by = join_by(patient_id))
  
  tmp_abx_days <- abx_days %>% 
    inner_join(tmp, by = join_by(patient_id)) %>% 
    filter(rx_date<index_date,
           rx_date>=(index_date-duration)) %>% 
    filter(abx_dates<index_date) %>% 
    distinct(patient_id,abx_dates) %>% 
    count(patient_id,name = "abx_days") 
  
  tmp_reg_data <- tmp %>% 
    select(patient_id,duration) %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_abx_days,by = join_by(patient_id)) %>% 
    mutate(abx_days=replace_na(abx_days,0L))  %>% 
    mutate(abx_days_cat = cut(abx_days,breaks = c(-1,0,3,5,7,10,99), labels = c("0","1-3","4-5","6-7","8-10",">10"))) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_days_cat, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}

rm(trial_num,tmp,out,tmp_reg_data,tmp_fit,tmp_abx_days)

mod_res3 <- distinct(sim_res,trial)

mod_res3 <- mod_res3 %>% 
  mutate(res = map(trial,run_dur_model3))

mod_res %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

mod_res2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

mod_res3 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))


######################################################
#### Models That Don't Restrict Based on Duration ####
######################################################

filter_abx_days <- abx_days %>% 
  inner_join(select(reg_data_base,patient_id,index_date)) %>% 
  filter(rx_date<index_date & rx_date>(index_date-29)) %>% 
  filter(abx_dates<index_date)

trial_num <- 1

run_dur_model1_2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num)
  
  tmp_inds <- filter_abx_days %>%
    distinct(patient_id) %>% 
    mutate(abx_vis=1L) %>% 
    inner_join(tmp,by = join_by(patient_id)) %>% 
    distinct(patient_id,abx_vis)
  
  tmp_reg_data <- tmp %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_inds,by = join_by(patient_id)) %>% 
    mutate(abx_vis=replace_na(abx_vis,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_vis, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}

mod_res1_2 <- distinct(sim_res,trial) %>% 
  mutate(res = map(trial,run_dur_model1_2))

mod_res1_2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

filter_abx_days_count <- filter_abx_days %>% 
  distinct(patient_id,abx_dates) %>% 
  count(patient_id)

run_dur_model2_2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num) 
  
  
  tmp_abx_days <- filter_abx_days_count %>% 
    inner_join(tmp, by = join_by(patient_id)) %>% 
    distinct(patient_id,abx_days=n) 
  
  tmp_reg_data <- tmp %>% 
    select(patient_id,duration) %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_abx_days,by = join_by(patient_id)) %>% 
    mutate(abx_days=replace_na(abx_days,0L)) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_days, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}

rm(trial_num,tmp,out,tmp_reg_data,tmp_fit,tmp_abx_days)

mod_res2_2 <- distinct(sim_res,trial) %>% 
  mutate(res = map(trial,run_dur_model2_2))

mod_res2_2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))

run_dur_model3_2 <- function(trial_num){
  
  tmp <- filter(sim_dur,trial==trial_num) 
  
  tmp_abx_days <- filter_abx_days_count %>% 
    inner_join(tmp, by = join_by(patient_id)) %>% 
    distinct(patient_id,abx_days=n) 
  
  tmp_reg_data <- tmp %>% 
    select(patient_id,duration) %>% 
    inner_join(reg_data_base,by = join_by(patient_id)) %>% 
    left_join(tmp_abx_days,by = join_by(patient_id)) %>% 
    mutate(abx_days=replace_na(abx_days,0L)) %>% 
    mutate(abx_days_cat = cut(abx_days,breaks = c(-1,0,3,5,7,10,99), labels = c("0","1-3","4-5","6-7","8-10",">10"))) 
  
  
  tmp_fit <- lm(duration~age_cat+sex+abx_days_cat, data = tmp_reg_data) 
  
  out <- coef(tmp_fit)
  
  return(out)
}

mod_res3_2 <- distinct(sim_res,trial) %>% 
  mutate(res = map(trial,run_dur_model3_2))

mod_res3_2 %>% 
  mutate(res = map(res,enframe)) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(mn_est = mean(value),
            lo_est = quantile(value,probs = c(0.025)),
            hi_est = quantile(value,probs = c(0.975)))



