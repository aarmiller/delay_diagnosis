
rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)
library(codeBuildr)

source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")

# name of condition
cond_name <- "sepsis_covid"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)

data_in_path <- paste0(delay_params$base_path,"delay_results/")
# sim_in_path <- paste0("/Volumes/AML/small_dbs/",cond_name,"/truven/enroll_restrict_365/","delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path,recursive = TRUE)
}

#### Load Data  ----------------------------------------------------------------

load(paste0(data_in_path,"all_dx_visits.RData"))
load(paste0(data_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes(delay_params$ssd_name) %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver) %>% 
  distinct() %>% 
  filter(!is.na(dx))

# load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

# extract patient ids and number of patients
patient_ids <- index_cases %>% 
  distinct(patient_id)

n_patients <- nrow(patient_ids)


# filter timemap down to final index cases
tm <- inner_join(tm, patient_ids, by = join_by(patient_id))

all_dx_visits <- inner_join(all_dx_visits, patient_ids, by = join_by("patient_id"))

ssd_dx_visits <- inner_join(all_dx_visits, ssd_codes, by = c("dx", "dx_ver"))


###############################
#### Compute Overall Trends####
###############################

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


####################
#### Fit Trends ####
####################

trends_out_path <- paste0(sim_out_path,"trends/")

if (!dir.exists(trends_out_path)) {
  dir.create(trends_out_path,recursive = TRUE)
}

### Fit all models

for (i in delay_params$cp){
  print(i)
  
  # fit all models -------------------------------------------------------------
  if (delay_params$periodicity){
    # if periodicity = TRUE add day of week term
    fit1_ssd <- lm(n ~ period + dow, filter(ssd_vis_count, period>i))
    fit2_ssd <- lm(n ~ log(period) + dow, filter(ssd_vis_count, period>i))
    fit3_ssd <- lm(n ~ poly(period,degree = 2) + dow, filter(ssd_vis_count, period>i))
    fit4_ssd <- lm(n ~ poly(period,degree = 3) + dow, filter(ssd_vis_count, period>i))
    
    fit1_all <- lm(n ~ period + dow, filter(all_vis_count, period>i))
    fit2_all <- lm(n ~ log(period) + dow, filter(all_vis_count, period>i))
    fit3_all <- lm(n ~ poly(period,degree = 2) + dow, filter(all_vis_count, period>i))
    fit4_all <- lm(n ~ poly(period,degree = 3) + dow, filter(all_vis_count, period>i))
  } else {
    fit1_ssd <- lm(n ~ period, filter(ssd_vis_count, period>i))
    fit2_ssd <- lm(n ~ log(period), filter(ssd_vis_count, period>i))
    fit3_ssd <- lm(n ~ poly(period,degree = 2), filter(ssd_vis_count, period>i))
    fit4_ssd <- lm(n ~ poly(period,degree = 3), filter(ssd_vis_count, period>i))
    
    fit1_all <- lm(n ~ period, filter(all_vis_count, period>i))
    fit2_all <- lm(n ~ log(period), filter(all_vis_count, period>i))
    fit3_all <- lm(n ~ poly(period,degree = 2), filter(all_vis_count, period>i))
    fit4_all <- lm(n ~ poly(period,degree = 3), filter(all_vis_count, period>i))
  }
  
  # compute rmse values --------------------------------------------------------
  mse_res_ssd <- tibble(model_name = c("fit1","fit2","fit3","fit4"),
                        model = c("linear","exponential","quadratic","cubic"),
                        rmse = c(sqrt(mean(fit1_ssd$residuals^2)),
                                 sqrt(mean(fit2_ssd$residuals^2)),
                                 sqrt(mean(fit3_ssd$residuals^2)),
                                 sqrt(mean(fit4_ssd$residuals^2)))) %>%
    mutate(label = paste("RMSE: ", round(rmse,2)))
  
  
  mse_res_all <- tibble(model_name = c("fit1","fit2","fit3","fit4"),
                        model = c("linear","exponential","quadratic","cubic"),
                        rmse = c(sqrt(mean(fit1_all$residuals^2)),
                                 sqrt(mean(fit2_all$residuals^2)),
                                 sqrt(mean(fit3_all$residuals^2)),
                                 sqrt(mean(fit4_all$residuals^2)))) %>%
    mutate(label = paste("RMSE: ", round(rmse,2)))
  
  # extract model fits ---------------------------------------------------------
  model_fits_ssd <- ssd_vis_count %>%
    mutate(linear = predict(fit1_ssd,newdata = .)) %>%
    mutate(exponential = predict(fit2_ssd,newdata = .)) %>%
    mutate(quadratic = predict(fit3_ssd,newdata = .)) %>%
    mutate(cubic = predict(fit4_ssd,newdata = .)) %>%
    gather(key = model, value = value, -period, -n,-dow) %>%
    inner_join(mse_res_ssd, by = "model") %>% 
    mutate(cp = i) %>% 
    mutate(num_miss = n-value) %>%
    mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
    mutate(num_miss = ifelse(period>cp,NA,num_miss)) 
  
  model_fits_all <- all_vis_count %>%
    mutate(linear = predict(fit1_all,newdata = .)) %>%
    mutate(exponential = predict(fit2_all,newdata = .)) %>%
    mutate(quadratic = predict(fit3_all,newdata = .)) %>%
    mutate(cubic = predict(fit4_all,newdata = .)) %>%
    gather(key = model, value = value, -period, -n,-dow) %>%
    inner_join(mse_res_all, by = "model") %>% 
    mutate(cp = i) %>% 
    mutate(num_miss = n-value) %>%
    mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
    mutate(num_miss = ifelse(period>cp,NA,num_miss)) 
  
  # Plot SSD results -----------------------------------------------------------
  # location to place RMSE values on plot
  y_pos <- .9*max(ssd_vis_count$n)
  x_pos <- .8*delay_params$upper_bound
  
  p1 <- model_fits_ssd %>%
    ggplot(aes(period,n)) +
    geom_line() +
    scale_x_reverse() +
    geom_line(aes(y = value), color = "red") +
    facet_wrap(~model) +
    geom_vline(aes(xintercept = i), linetype =2) +
    theme_bw() +
    geom_text(data = mse_res_ssd,
              mapping = aes(x = x_pos, y = y_pos, label = label))
  
  ggsave(filename = paste0(trends_out_path,"expected_trends_ssd_cp",i,".pdf"), plot = p1)
  
  p2 <- model_fits_ssd %>%
    filter(period<=cp) %>%
    ggplot(aes(period,num_miss)) +
    geom_histogram(stat = "identity") +
    scale_x_reverse() +
    ylab("Number of Missed Opportunities") +
    xlab("Days Before Index Diagnosis") +
    theme_bw() +
    facet_wrap(~model) +
    ggtitle(paste0("Estimated Number of Missed Opportunities for ",cond_name))
  
  ggsave(filename =  paste0(trends_out_path,"num_misses_ssd_cp",i,".pdf"), plot = p2)
  
  # Plat all visit results -----------------------------------------------------
  # location to place RMSE values on plot
  y_pos <- .9*max(all_vis_count$n)
  x_pos <- .8*delay_params$upper_bound
  
  p1 <- model_fits_all %>%
    ggplot(aes(period,n)) +
    geom_line() +
    scale_x_reverse() +
    geom_line(aes(y = value), color = "red") +
    facet_wrap(~model) +
    geom_vline(aes(xintercept = i), linetype =2) +
    theme_bw() +
    geom_text(data = mse_res_all,
              mapping = aes(x = x_pos, y = y_pos, label = label))
  
  ggsave(filename = paste0(trends_out_path,"expected_trends_all_cp",i,".pdf"), plot = p1)
  
  p2 <- model_fits_all %>%
    filter(period<=cp) %>%
    ggplot(aes(period,num_miss)) +
    geom_histogram(stat = "identity") +
    scale_x_reverse() +
    ylab("Number of Missed Opportunities") +
    xlab("Days Before Index Diagnosis") +
    theme_bw() +
    facet_wrap(~model) +
    ggtitle(paste0("Estimated Number of Missed Opportunities for ",cond_name))
  
  ggsave(filename =  paste0(trends_out_path,"num_misses_all_cp",i,".pdf"), plot = p2)
  
}
rm(fit1_all,fit2_all,fit3_all,fit4_all,
   fit1_ssd,fit2_ssd,fit3_ssd,fit4_ssd,
   model_fits_all,model_fits_ssd,p1,p2,
   mse_res_ssd,mse_res_all,y_pos,x_pos,i)

#### fit final trends ----------------------------------------------------------

models <- tibble(cp = delay_params$cp) %>% 
  mutate(model = map(cp,~delay_params$final_model)) %>% 
  unnest(model)


# all visits
all_vis_count_fitted <- models %>% 
  mutate(counts = map2(model,cp,
                       ~return_fits(data = all_vis_count,
                                    model = .x,
                                    cp = .y,
                                    periodicity = delay_params$periodicity))) %>% 
  mutate(num_miss=map(counts,~summarise(.,num_miss = sum(num_miss,na.rm = TRUE)))) %>% 
  unnest(num_miss)

# ssd visits
ssd_vis_count_fitted <- models %>% 
  mutate(counts = map2(model,cp,
                       ~return_fits(data = ssd_vis_count,
                                    model = .x,
                                    cp = .y,
                                    periodicity = delay_params$periodicity))) %>% 
  mutate(num_miss=map(counts,~summarise(.,num_miss = sum(num_miss,na.rm = TRUE)))) %>% 
  unnest(num_miss)

save(all_vis_count_fitted, ssd_vis_count_fitted,
     file = paste0(trends_out_path,"fit_trends.RData"))



################################
#### Prepare Bootstrap Data ####
################################

#### Setup outer bootstrap data ------------------------------------------------

# sample overall set of patients for each bootstrap
set.seed(192837)
tmp <- tibble(boot_trial=1:delay_params$boot_trials) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(patient_id = sample(patient_ids$patient_id, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# tmp$boot_sample[[1]]

#### compute counts for all visits ---------------------------------------------
get_all_counts <- function(patient_set){
  patient_set %>% 
    inner_join(all_dx_visits, by = join_by("patient_id"),relationship = "many-to-many") %>%
    mutate(period = -days_since_index) %>%
    distinct(patient_id,boot_id,period) %>% 
    count(period) %>% 
    filter(period>0) %>% 
    mutate(dow = as.factor(period %% 7))
}

# get_all_counts(tmp$boot_sample[[1]])


# For each set of patients selected pull the time map visit counts
tmp <- tmp %>% 
  mutate(all_vis_count = map(boot_sample,get_all_counts))



#### compute counts for ssd visits ---------------------------------------------
get_ssd_counts <- function(patient_set){
  patient_set %>% 
    inner_join(ssd_dx_visits, by = join_by("patient_id"),relationship = "many-to-many") %>%
    mutate(period = -days_since_index) %>%
    distinct(patient_id,boot_id,period) %>% 
    count(period) %>% 
    filter(period>0) %>% 
    mutate(dow = as.factor(period %% 7))
}

# get_ssd_counts(tmp$boot_sample[[1]])

# for each set of patients pull the ssd counts
tmp <- tmp %>% 
  mutate(ssd_vis_count = map(boot_sample,get_ssd_counts))

# tmp$all_vis_count[[1]]
# tmp$ssd_vis_count[[1]]


### fit models -----------------------------------------------------------------
fit_models <- function(count_data){
  models %>% 
    mutate(counts = map2(model,cp,
                         ~return_fits(data = count_data,
                                      model = .x,
                                      cp = .y,
                                      periodicity = delay_params$periodicity)))
}

# fit trends for both all visits and ssd visits
tmp <- tmp %>% 
  mutate(all_vis_count = map(all_vis_count,fit_models)) %>% 
  mutate(ssd_vis_count = map(ssd_vis_count,fit_models))

boot_data <- tmp

# compute number missed
tmp1 <- boot_data %>% 
  select(-boot_sample) %>% 
  unnest(all_vis_count) %>% 
  mutate(n_miss_all = map(counts,~summarise(.,n_miss_all = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_all) %>% 
  select(boot_trial,cp,model,n_miss_all)

tmp2 <- boot_data %>% 
  select(-boot_sample) %>% 
  unnest(ssd_vis_count) %>% 
  mutate(n_miss_ssd = map(counts,~summarise(.,n_miss_ssd = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_ssd) %>% 
  select(boot_trial,cp,model,n_miss_ssd)

boot_miss_counts <- inner_join(tmp1,tmp2)

#### finalize and save data for bootstrapping ----------------------------------
boot_data <- tmp
rm(tmp,tmp1,tmp2)

# save bootstrap data
save(boot_data,boot_miss_counts,file = paste0(sim_out_path,"boot_data.RData"))


###########################################
#### Run bootstrap simulation analysis ####
###########################################

cluster <- parallel::makeCluster(30)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterCall(cluster, function() library(delaySim))
parallel::clusterExport(cluster,c("run_boot_trials","boot_data","sim_tm_all","sim_tm_ssd","delay_params","n_patients"),
                        envir=environment())


for (i in 1:nrow(models)) {
  use_model <- models$model[i]
  use_cp <- models$cp[i]
  
  parallel::clusterExport(cluster,c("use_model","use_cp"))
  
  print(paste0("Model: ",use_model,"; CP: ",use_cp))
  
  out_path <- paste0(sim_out_path,use_model,"_cp",use_cp,"/")
  # print(out_path)
  
  if (!dir.exists(out_path)) { dir.create(out_path) }
  
  set.seed(5678)
  sim_res_all <- parallel::parLapply(cl = cluster,
                                     1:delay_params$boot_trials,
                                     function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                                 tm_data = sim_tm_all,
                                                                 miss_bins = filter(boot_data$all_vis_count[[x]],
                                                                                    cp == use_cp, model == use_model) %>%
                                                                   .$counts %>%
                                                                   .[[1]] %>% 
                                                                   select(period,num_miss),
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
                                                                 miss_bins = filter(boot_data$ssd_vis_count[[x]],
                                                                                    cp == use_cp, model == use_model) %>%
                                                                   .$counts %>%
                                                                   .[[1]] %>%
                                                                   select(period,num_miss),
                                                                 delay_params = delay_params,
                                                                 n_patients = n_patients,
                                                                 n_trials = delay_params$sim_trials)})

  sim_res_ssd <- map2(sim_res_ssd,1:delay_params$boot_trials,~mutate(.x,boot_trial=.y)) %>%
    bind_rows()

  save(sim_res_all,file = paste0(out_path,"sim_res_all.RData"))
  save(sim_res_ssd,file = paste0(out_path,"sim_res_ssd.RData"))
  
  rm(sim_res_all,sim_res_ssd)
  
}


#### Extract sim result observations -------------------------------------------

# # now extract the specific sim obs values for future steps
# tmp1 <- sim_res_ssd %>% unnest(res) %>% distinct(obs)
# tmp2 <- sim_res_all %>% unnest(res) %>% distinct(obs)
# 
# sim_obs_reduced <- sim_obs %>% 
#   inner_join(distinct(bind_rows(tmp1,tmp2)), by = "obs") %>% 
#   mutate(period = -days_since_index) %>% 
#   select(obs,patient_id,period)
# 
# save(sim_obs_reduced, file = paste0(sim_out_path,"sim_obs_reduced.RData"))
# 
# # cleanup
# rm(list = ls()[!(ls() %in% c("delay_params","sim_out_path","sim_in_path","cond_name","n_patients"))])
# gc()


