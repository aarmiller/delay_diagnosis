
rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(furrr)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "dengue"

# load parameters
load("/Shared/AML/params/delay_any_params.RData")
delay_params <- delay_any_params[[cond_name]]

# input path
sim_in_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")

# output path
sim_out_path <- paste0(sim_in_path,"sim_results/")

# create output directory
if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

# number of trials and bootstraps
n_trials <- 10
n_bootstraps <- 100

##############################
#### Simulation Functions ####
##############################
source("github/delay_diagnosis/build_scripts/R/functions/simulation_functions.R")
source("github/delay_diagnosis/build_scripts/R/functions/trend_functions.R")
# source("build_scripts/R/functions/simulation_functions.R")


######################
#### Prepare Data ####
######################

#### Load Data  ----------------------------------------------------------------

load(paste0(sim_in_path,"all_dx_visits.RData"))
load(paste0(sim_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver) %>% 
  distinct()

# index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

sim_tm_ssd <- all_dx_visits %>%
  mutate(period = -days_since_index) %>%
  inner_join(ssd_codes,by = c("dx", "dx_ver")) %>%
  distinct(patient_id,period,days_since_index) %>%
  inner_join(sim_obs,by = c("patient_id", "days_since_index")) 

sim_tm_any <- all_dx_visits %>%
  mutate(period = -days_since_index) %>%
  distinct(patient_id,period,days_since_index) %>%
  left_join(sim_obs,.,by = c("patient_id", "days_since_index"))

all_vis_count_ssd <- sim_tm_ssd %>%
  count(period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

all_vis_count_any <- sim_tm_any %>%
  count(period) %>%
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))


# # merge observations into time map to extract visit types
obs_tm_ssd <- sim_tm_ssd %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,obs_stay,inpatient)

obs_tm_any <- sim_tm_any %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,obs_stay,inpatient)

# extract patient ids and number of patients
patient_ids <- index_dx_dates %>% distinct(patient_id)
n_patients <- nrow(patient_ids)

#############################
#### Fit Expected Trends ####
#############################

#### fit SSD models ------------------------------------------------------------
if (delay_params$periodicity){
  # if periodicity = TRUE add day of week term
  fit1 <- lm(n ~ period + dow, filter(all_vis_count_ssd, period>=delay_params$cp))
  fit2 <- lm(n ~ log(period) + dow, filter(all_vis_count_ssd, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2) + dow, filter(all_vis_count_ssd, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3) + dow, filter(all_vis_count_ssd, period>delay_params$cp))
} else {
  fit1 <- lm(n ~ period, filter(all_vis_count_ssd, period>delay_params$cp))
  fit2 <- lm(n ~ log(period), filter(all_vis_count_ssd, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2), filter(all_vis_count_ssd, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3), filter(all_vis_count_ssd, period>delay_params$cp))
}

# compute rmse values
mse_res_ssd <- tibble(model_name = c("fit1","fit2","fit3","fit4"),
                  model = c("linear","exponential","quadratic","cubic"),
                  rmse = c(sqrt(sum(fit1$residuals^2)),
                           sqrt(sum(fit2$residuals^2)),
                           sqrt(sum(fit3$residuals^2)),
                           sqrt(sum(fit4$residuals^2)))) %>%
  mutate(label = paste("RMSE: ", round(rmse,2)))

# location to place RMSE values on plot
y_pos <- .9*max(all_vis_count_ssd$n)
x_pos <- .8*delay_params$upper_bound

model_fits_ssd <- all_vis_count_ssd %>%
  mutate(linear = predict(fit1,newdata = .)) %>%
  mutate(exponential = predict(fit2,newdata = .)) %>%
  mutate(quadratic = predict(fit3,newdata = .)) %>%
  mutate(cubic = predict(fit4,newdata = .)) %>%
  gather(key = model, value = value, -period, -n,-dow) %>%
  inner_join(mse_res_ssd, by = "model") %>% 
  mutate(cp = delay_params$cp) %>% 
  mutate(num_miss = n-value) %>%
  mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
  mutate(num_miss = ifelse(period>cp,NA,num_miss)) 

## plot trends
p1 <- model_fits_ssd %>%
  ggplot(aes(period,n)) +
  geom_line() +
  scale_x_reverse() +
  geom_line(aes(y = value), color = "red") +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = delay_params$cp), linetype =2) +
  theme_bw() +
  geom_text(data = mse_res_ssd,
            mapping = aes(x = x_pos, y = y_pos, label = label))

ggsave(filename = paste0(sim_out_path,"/expected_trends_ssd_visits.pdf"),
       plot = p1)

## Plot number of missed opportunities
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

ggsave(filename = paste0(sim_out_path,"/number_of_missed_opportunities_ssd_visits.pdf"),
       plot = p2)

#### fit Any Visit models ------------------------------------------------------
if (delay_params$periodicity){
  # if periodicity = TRUE add day of week term
  fit1 <- lm(n ~ period + dow, filter(all_vis_count_any, period>=delay_params$cp))
  fit2 <- lm(n ~ log(period) + dow, filter(all_vis_count_any, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2) + dow, filter(all_vis_count_any, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3) + dow, filter(all_vis_count_any, period>delay_params$cp))
} else {
  fit1 <- lm(n ~ period, filter(all_vis_count_any, period>delay_params$cp))
  fit2 <- lm(n ~ log(period), filter(all_vis_count_any, period>delay_params$cp))
  fit3 <- lm(n ~ poly(period,degree = 2), filter(all_vis_count_any, period>delay_params$cp))
  fit4 <- lm(n ~ poly(period,degree = 3), filter(all_vis_count_any, period>delay_params$cp))
}

# compute rmse values
mse_res_any <- tibble(model_name = c("fit1","fit2","fit3","fit4"),
                      model = c("linear","exponential","quadratic","cubic"),
                      rmse = c(sqrt(sum(fit1$residuals^2)),
                               sqrt(sum(fit2$residuals^2)),
                               sqrt(sum(fit3$residuals^2)),
                               sqrt(sum(fit4$residuals^2)))) %>%
  mutate(label = paste("RMSE: ", round(rmse,2)))

# location to place RMSE values on plot
y_pos <- .9*max(all_vis_count_any$n)
x_pos <- .8*delay_params$upper_bound

model_fits_any <- all_vis_count_any %>%
  mutate(linear = predict(fit1,newdata = .)) %>%
  mutate(exponential = predict(fit2,newdata = .)) %>%
  mutate(quadratic = predict(fit3,newdata = .)) %>%
  mutate(cubic = predict(fit4,newdata = .)) %>%
  gather(key = model, value = value, -period, -n,-dow) %>%
  inner_join(mse_res_any, by = "model") %>% 
  mutate(cp = delay_params$cp) %>% 
  mutate(num_miss = n-value) %>%
  mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
  mutate(num_miss = ifelse(period>cp,NA,num_miss)) 

## plot trends
p1 <- model_fits_any %>%
  ggplot(aes(period,n)) +
  geom_line() +
  scale_x_reverse() +
  geom_line(aes(y = value), color = "red") +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = delay_params$cp), linetype =2) +
  theme_bw() +
  geom_text(data = mse_res_any,
            mapping = aes(x = x_pos, y = y_pos, label = label))

ggsave(filename = paste0(sim_out_path,"/expected_trends_any_visit.pdf"),
       plot = p1)

## Plot number of missed opportunities
p2 <- model_fits_any %>%
  filter(period<=cp) %>%
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Diagnosis") +
  theme_bw() +
  facet_wrap(~model) +
  ggtitle(paste0("Estimated Number of Missed Opportunities for ",cond_name))

ggsave(filename = paste0(sim_out_path,"/number_of_missed_opportunities_any_visit.pdf"),
       plot = p2)

################################
#### Prepare Bootstrap Data ####
################################

#### Extract Optimal Models ----------------------------------------------------

if (is.na(delay_params$final_model)) {
  selected_model_ssd <- mse_res_ssd %>%
    filter(model!="exponential") %>%
    filter(rmse ==min(rmse))
  
  selected_model_any <- mse_res_any %>%
    filter(model!="exponential") %>%
    filter(rmse ==min(rmse))
  
} else {
  selected_model_ssd <- mse_res_ssd %>%
    filter(model==delay_params$final_model) %>% 
    filter(rmse ==min(rmse))
  
  selected_model_any <- mse_res_any %>%
    filter(model==delay_params$final_model) %>% 
    filter(rmse ==min(rmse))
}

#### Setup outer bootstrap data ------------------------------------------------

# sample overall set of patients for each bootstrap
tmp <- tibble(boot_trial=1:n_bootstraps) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(patient_id = sample(patient_ids$patient_id, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# tmp$boot_sample[[1]]

#### compute counts for all visits ---------------------------------------------
# For each set of patients selected pull the time map visit counts
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "patient_id",relationship = "many-to-many") %>%
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
                      ~inner_join(.,all_dx_visits, by = "patient_id",relationship = "many-to-many") %>%
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
                                                        model = selected_model_any$model,
                                                        cp = delay_params$cp,
                                                        periodicity = delay_params$periodicity)),
         ssd_vis_count = map(ssd_vis_count,~return_fits(.,
                                                        model = selected_model_ssd$model,
                                                        cp = delay_params$cp,
                                                        periodicity = delay_params$periodicity)))
# compute number missed
tmp <- tmp %>% 
  mutate(n_miss_all = map(all_vis_count,
                          ~filter(., period<delay_params$cp) %>%  
                            summarise(.,n_miss_all = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_all) %>% 
  mutate(n_miss_ssd = map(ssd_vis_count,
                          ~filter(., period<delay_params$cp) %>%  
                            summarise(.,n_miss_ssd = sum(round(num_miss,0),na.rm = T)))) %>% 
  unnest(n_miss_ssd) %>% 
  mutate(n_expected_all =  map(all_vis_count,
                               ~filter(., period<delay_params$cp) %>%  
                                 summarise(.,n_expected_all = sum(round(pred,0),na.rm = T)))) %>% 
  unnest(n_expected_all) %>% 
  mutate(n_expected_ssd =  map(ssd_vis_count,
                               ~filter(., period<delay_params$cp) %>%  
                                 summarise(.,n_expected_ssd = sum(round(pred,0),na.rm = T)))) %>% 
  unnest(n_expected_ssd)  %>% 
  mutate(n_total_all =  map(all_vis_count,
                            ~filter(., period<delay_params$cp) %>%  
                              summarise(.,n_total_all = sum(round(n,0),na.rm = T)))) %>% 
  unnest(n_total_all) %>% 
  mutate(n_total_ssd =  map(ssd_vis_count,
                            ~filter(., period<delay_params$cp) %>%  
                              summarise(.,n_total_ssd = sum(round(n,0),na.rm = T)))) %>% 
  unnest(n_total_ssd) 

boot_data <- tmp

rm(tmp)

save(boot_data, file = paste0(sim_out_path,"boot_data.RData"))

###########################################
#### Run bootstrap simulation analysis ####
###########################################

cluster <- parallel::makeCluster(30)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterCall(cluster, function() library(delaySim))
parallel::clusterExport(cluster,c("run_boot_trials","boot_data","sim_tm_any","sim_tm_ssd","delay_params","n_patients","n_trials"),
                        envir=environment())

set.seed(5678)
sim_res_any <- parallel::parLapply(cl = cluster,
                                   1:n_bootstraps,
                                   function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                               tm_data = sim_tm_any,
                                                               miss_bins = select(boot_data$all_vis_count[[x]],period,num_miss),
                                                               delay_params = delay_params,
                                                               n_patients = n_patients,
                                                               n_trials = n_trials)})

sim_res_any <- map2(sim_res_any,1:n_bootstraps,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

sim_res_any$res[1]

set.seed(5678)
sim_res_ssd <- parallel::parLapply(cl = cluster,
                                   1:n_bootstraps,
                                   function(x){run_boot_trials(boot_ids = boot_data$boot_sample[[x]],
                                                               tm_data = sim_tm_ssd,
                                                               miss_bins = select(boot_data$ssd_vis_count[[x]],period,num_miss),
                                                               delay_params = delay_params,
                                                               n_patients = n_patients,
                                                               n_trials = n_trials)})

sim_res_ssd <- map2(sim_res_ssd,1:n_bootstraps,~mutate(.x,boot_trial=.y)) %>% 
  bind_rows()

parallel::stopCluster(cluster)

rm(cluster)

## Save Simulation Results -----------------------------------------------------

save(sim_res_ssd, sim_tm_ssd, obs_tm_ssd, file = paste0(sim_out_path,"sim_res_ssd.RData"))
save(sim_res_any, sim_tm_any, obs_tm_any, file = paste0(sim_out_path,"sim_res_any.RData"))




#################################
#### Generate Setting Counts ####
#################################

load(paste0(sim_in_path,"caseids.RData"))


## Compute Index Setting Counts ------------------------------------------------

out_index_n <- tm %>% 
  filter(days_since_index==0) %>% 
  distinct(patient_id, outpatient,ed,inpatient,obs_stay) %>% 
  gather(key = Setting, value = value, -patient_id) %>% 
  group_by(Setting) %>% 
  summarise(index_n = sum(value))

out_index_n <- bind_rows(out_index_n,
                         filter(out_index_n,Setting=="inpatient") %>% 
                           mutate(Setting = "inpatient visit"))

## Compute Potential Opportunity Counts - SSD ----------------------------------

# Count for visit days
out_pot_opps_ssd <- all_dx_visits %>% 
  inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  filter(between(days_since_index,-delay_params$cp+1,-1)) %>% ### NEED TO DOUBLE CHECK THIS WITH THE NEW ALIGNMENT                  
  distinct(obs,outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value, -obs) %>%
  group_by(Setting) %>% 
  summarise(potential_opps = sum(value))

# build counts for inpatient stays
tmp <- all_dx_visits %>% 
  inner_join(ssd_codes,by = join_by(dx, dx_ver)) %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
  inner_join(distinct(caseids,obs,caseid,patient_id),by = join_by(patient_id, obs)) %>% 
  filter(between(days_since_index,-delay_params$cp+1,-1)) %>% 
  distinct(patient_id,caseid) %>% 
  count(name = "potential_opps") %>% 
  mutate(Setting = "inpatient visit")

# combine daily settings and inpatient visits
out_pot_opps_ssd <- bind_rows(out_pot_opps_ssd,tmp)

## Compute Potential Opportunity Counts - Any ----------------------------------

# Count for visit days
out_pot_opps_any <- all_dx_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  filter(between(days_since_index,-delay_params$cp+1,-1)) %>% ### NEED TO DOUBLE CHECK THIS WITH THE NEW ALIGNMENT                  
  distinct(obs,outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value, -obs) %>%
  group_by(Setting) %>% 
  summarise(potential_opps = sum(value))

# build counts for inpatient stays
tmp <- all_dx_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
  inner_join(distinct(caseids,obs,caseid,patient_id),by = join_by(patient_id, obs)) %>% 
  filter(between(days_since_index,-delay_params$cp+1,-1)) %>% 
  distinct(patient_id,caseid) %>% 
  count(name = "potential_opps") %>% 
  mutate(Setting = "inpatient visit")

# combine daily settings and inpatient visits
out_pot_opps_any <- bind_rows(out_pot_opps_any,tmp)

## Compute Setting Miss Counts - SSD -------------------------------------------

obs_tm <- sim_obs %>%
  distinct(obs,days_since_index,patient_id) %>%
  inner_join(tm,by = c("days_since_index", "patient_id")) %>%
  distinct(obs,outpatient,ed,obs_stay,inpatient)

# # merge observations into time map to extract visit types
tmp <- sim_res_ssd %>% 
  mutate(res = map(res,~inner_join(.,obs_tm,by = "obs")))

# get setting counts for each trial
setting_counts_ssd <- tmp %>% 
  mutate(n = map_int(res,nrow)) %>% 
  mutate(outpatient = map_int(res,~sum(.$outpatient))) %>% 
  mutate(ed = map_int(res,~sum(.$ed))) %>% 
  mutate(obs_stay = map_int(res,~sum(.$obs_stay))) %>% 
  mutate(inpatient = map_int(res,~sum(.$inpatient))) %>% 
  select(sim_trial,boot_trial,n:inpatient)

# compute miss counts for outpatient, ed, inpatient days and obs_stay
out1 <- setting_counts_ssd %>% 
  select(-n) %>% 
  gather(key = Setting, value = n, -sim_trial,-boot_trial) %>% 
  group_by(Setting) %>% 
  summarise(miss_mean = round(mean(n),0),
            miss_median = round(median(n),0),
            miss_low = round(quantile(n,probs = 0.025),0),
            miss_high = round(quantile(n,probs = 0.975),0))

# merge caseids into simulation results
tmp <- sim_res_ssd %>% 
  mutate(res = map(res,~inner_join(.,distinct(caseids,obs,caseid,patient_id),by = "obs",relationship = "many-to-many") %>% 
                     distinct(patient_id,boot_id,caseid)))

# compute counts for inpatient stays
out2 <- tmp %>%
  mutate(n_vis = map_int(res,nrow)) %>% 
  summarise(miss_mean = round(mean(n_vis),0),
            miss_median = round(median(n_vis),0),
            miss_low = round(quantile(n_vis,probs = 0.025),0),
            miss_high = round(quantile(n_vis,probs = 0.975),0)) %>% 
  mutate(`Setting`= "inpatient visit")

out_miss_n_ssd <- bind_rows(out1,out2)

## Compute Setting Miss Counts - Any -------------------------------------------

# # merge observations into time map to extract visit types
tmp <- sim_res_any %>% 
  mutate(res = map(res,~inner_join(.,obs_tm,by = "obs")))

# get setting counts for each trial
setting_counts_any <- tmp %>% 
  mutate(n = map_int(res,nrow)) %>% 
  mutate(outpatient = map_int(res,~sum(.$outpatient))) %>% 
  mutate(ed = map_int(res,~sum(.$ed))) %>% 
  mutate(obs_stay = map_int(res,~sum(.$obs_stay))) %>% 
  mutate(inpatient = map_int(res,~sum(.$inpatient))) %>% 
  select(sim_trial,boot_trial,n:inpatient)

# compute miss counts for outpatient, ed, inpatient days and obs_stay
out1 <- setting_counts_any %>% 
  select(-n) %>% 
  gather(key = Setting, value = n, -sim_trial,-boot_trial) %>% 
  group_by(Setting) %>% 
  summarise(miss_mean = round(mean(n),0),
            miss_median = round(median(n),0),
            miss_low = round(quantile(n,probs = 0.025),0),
            miss_high = round(quantile(n,probs = 0.975),0))

# merge caseids into simulation results
tmp <- sim_res_any %>% 
  mutate(res = map(res,~inner_join(.,distinct(caseids,obs,caseid,patient_id),by = "obs",relationship = "many-to-many") %>% 
                     distinct(patient_id,boot_id,caseid)))

# compute counts for inpatient stays
out2 <- tmp %>%
  mutate(n_vis = map_int(res,nrow)) %>% 
  summarise(miss_mean = round(mean(n_vis),0),
            miss_median = round(median(n_vis),0),
            miss_low = round(quantile(n_vis,probs = 0.025),0),
            miss_high = round(quantile(n_vis,probs = 0.975),0)) %>% 
  mutate(`Setting`= "inpatient visit")

out_miss_n_any <- bind_rows(out1,out2)

out_miss_n_any
out_miss_n_ssd

### Compute Percentage Miss ----------------------------------------------------

setting_counts_any

out_miss_frac_ssd <- setting_counts_ssd %>% 
  mutate_at(vars(outpatient,ed,obs_stay,inpatient),~round(100*./n,2)) %>% 
  select(outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean = mean(value),
            frac_median = median(value),
            frac_low = round(quantile(value,probs = 0.025),2),
            frac_high = round(quantile(value,probs = 0.975),2))

out_miss_frac_any <- setting_counts_any %>% 
  mutate_at(vars(outpatient,ed,obs_stay,inpatient),~round(100*./n,2)) %>% 
  select(outpatient,ed,obs_stay,inpatient) %>% 
  gather(key = Setting, value = value) %>% 
  group_by(Setting) %>% 
  summarise(frac_mean = mean(value),
            frac_median = median(value),
            frac_low = round(quantile(value,probs = 0.025),2),
            frac_high = round(quantile(value,probs = 0.975),2))

## Assemble final setting tables -----------------------------------------------

# Final Setting Count Table
setting_table_ssd <- tibble(Setting = c("outpatient","ed","obs_stay","inpatient","inpatient visit")) %>% 
  left_join(out_index_n, by = "Setting") %>% 
  left_join(out_pot_opps_ssd, by = "Setting") %>% 
  left_join(out_miss_n_ssd, by = "Setting") %>% 
  left_join(out_miss_frac_ssd, by = "Setting")

setting_table_any <- tibble(Setting = c("outpatient","ed","obs_stay","inpatient","inpatient visit")) %>% 
  left_join(out_index_n, by = "Setting") %>% 
  left_join(out_pot_opps_any, by = "Setting") %>% 
  left_join(out_miss_n_any, by = "Setting") %>% 
  left_join(out_miss_frac_any, by = "Setting")

# Output Table for Reporting
# setting_table_any %>% 
#   mutate(`Miss Counts` = paste0(miss_mean, " (", miss_low, "-", miss_high, ")")) %>% 
#   mutate(`Miss Percentage` = paste0(round(frac_mean,2), " (", frac_low, "-", frac_high, ")")) %>% 
#   select(Setting,`Index Count`=index_n, `Potential Opportunities`=potential_opps,
#          `Miss Counts`, `Miss Percentage`) %>% 
#   mutate(`Miss Percentage` = ifelse(Setting == "inpatient visit", "", `Miss Percentage`))
# 
# setting_table_ssd %>% 
#   mutate(`Miss Counts` = paste0(miss_mean, " (", miss_low, "-", miss_high, ")")) %>% 
#   mutate(`Miss Percentage` = paste0(round(frac_mean,2), " (", frac_low, "-", frac_high, ")")) %>% 
#   select(Setting,`Index Count`=index_n, `Potential Opportunities`=potential_opps,
#          `Miss Counts`, `Miss Percentage`) %>% 
#   mutate(`Miss Percentage` = ifelse(Setting == "inpatient visit", "", `Miss Percentage`))


##########################################
#### Compute Results for Delay Report ####
##########################################

sim_results_output <- list()

## General Information ---------------------------------------------------------

# number of patients
sim_results_output$n_patients <- n_patients

# number of trials and bootstraps
sim_results_output$n_trials <- n_trials
sim_results_output$n_bootstraps <- n_bootstraps

# change-point and optimal model
sim_results_output$cp <- delay_params$cp
sim_results_output$selected_model_ssd <- selected_model_ssd
sim_results_output$selected_model_any <- selected_model_any

# aggregate potential opportunities
sim_results_output$pot_opps_ssd <- nrow(filter(sim_tm_ssd, between(period,1,delay_params$cp)))
sim_results_output$pot_opps_any <- nrow(filter(sim_tm_any, between(period,1,delay_params$cp)))

## Results for generating plots ------------------------------------------------

# overall model fits
sim_results_output$model_fits_ssd <- model_fits_ssd
sim_results_output$model_fits_any <- model_fits_any

sim_results_output$mse_res_ssd <- mse_res_ssd
sim_results_output$mse_res_any <- mse_res_any

# bootstrap model fits
sim_results_output$boot_fits <- select(boot_data,-boot_sample)

## Compute Simulation Statistics -----------------------------------------------

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

plan(multisession, workers = 30)

# SSD results
tmp1 <- sim_res_ssd %>%
  mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
                                                         sim_tm_ssd,
                                                         delay_params = delay_params,
                                                         n_patients = n_patients)))

tmp1 <- tmp1 %>%
  select(sim_trial,boot_trial,sim_stats) %>%
  mutate(stats = map(sim_stats,~.$stats)) %>%
  mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>%
  mutate(dur_tab = map(sim_stats,~.$dur_tab))

sim_stats_ssd <- tmp1 %>%
  select(-sim_stats)

# Any visit results
tmp1 <- sim_res_any %>%
  mutate(sim_stats = future_map(res, ~compute_boot_stats(.,
                                                         sim_tm_any,
                                                         delay_params = delay_params,
                                                         n_patients = n_patients)))

tmp1 <- tmp1 %>%
  select(sim_trial,boot_trial,sim_stats) %>%
  mutate(stats = map(sim_stats,~.$stats)) %>%
  mutate(miss_tab = map(sim_stats,~.$miss_tab)) %>%
  mutate(dur_tab = map(sim_stats,~.$dur_tab))

sim_stats_any <- tmp1 %>%
  select(-sim_stats)

# save base sim stats
save(sim_stats_ssd,sim_stats_any,file = paste0(sim_out_path,"sim_stats.RData"))
rm(tmp1)

## Aggregate Simulation Statistics ---------------------------------------------

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

agg_stats_ssd <- aggregate_sim_stats(sim_stats_ssd)
agg_stats_any <- aggregate_sim_stats(sim_stats_any)

aggregate_miss_tab <- function(sim_stats_data){
  tmp1 <- sim_stats_data %>%
    select(sim_trial,boot_trial,miss_tab) %>%
    unnest(miss_tab) %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~round(mean(.),2)) %>%
    gather(key = key, value = mean, -miss_bin) %>%
    mutate(mean = round(mean,2))
  
  tmp2 <- sim_stats_data %>%
    select(sim_trial,boot_trial,miss_tab) %>%
    unnest(miss_tab) %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>%
    gather(key = key, value = low, -miss_bin) %>%
    mutate(low = round(low,2))
  
  tmp3 <- sim_stats_data %>%
    select(sim_trial,boot_trial,miss_tab) %>%
    unnest(miss_tab) %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>%
    gather(key = key, value = high, -miss_bin) %>%
    mutate(high = round(high,2))
  
  miss_bins_output <- inner_join(tmp1,tmp2) %>%
    inner_join(tmp3) %>%
    mutate(out = paste0(mean," (",low,"-",high,")")) %>%
    select(miss_bin,key,out) %>%
    spread(key = key, value = out)
  
  return(miss_bins_output)
}

agg_miss_bins_ssd <- aggregate_miss_tab(sim_stats_ssd)
agg_miss_bins_any <- aggregate_miss_tab(sim_stats_any)

aggregate_duration_tab <- function(sim_stats_data){
  tmp1 <- sim_stats_data %>%
    select(sim_trial,boot_trial,dur_tab) %>%
    unnest(dur_tab) %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~round(mean(.),2)) %>%
    gather(key = key, value = mean, -duration_bin) %>%
    mutate(mean = round(mean,2))
  
  tmp2 <- sim_stats_data %>%
    select(sim_trial,boot_trial,dur_tab) %>%
    unnest(dur_tab) %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.025)),2)) %>%
    gather(key = key, value = low, -duration_bin) %>%
    mutate(low = round(low,2))
  
  tmp3 <- sim_stats_data %>%
    select(sim_trial,boot_trial,dur_tab) %>%
    unnest(dur_tab) %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~round(quantile(.,probs = c(0.975)),2)) %>%
    gather(key = key, value = high, -duration_bin) %>%
    mutate(high = round(high,2))
  
  dur_bins_output <- inner_join(tmp1,tmp2) %>%
    inner_join(tmp3) %>%
    mutate(out = paste0(mean," (",low,"-",high,")")) %>%
    select(duration_bin,key,out) %>%
    spread(key = key, value = out)
  
  return(dur_bins_output)
}


agg_duration_bins_ssd <- aggregate_duration_tab(sim_stats_ssd)
agg_duration_bins_any <- aggregate_duration_tab(sim_stats_any)

# add to output
sim_results_output$agg_stats_ssd <- agg_stats_ssd 
sim_results_output$agg_stats_any <- agg_stats_any
sim_results_output$agg_miss_bins_ssd <- agg_miss_bins_ssd
sim_results_output$agg_miss_bins_any <- agg_miss_bins_any
sim_results_output$agg_duration_bins_ssd <- agg_duration_bins_ssd
sim_results_output$agg_duration_bins_any <- agg_duration_bins_any

## Add setting information -----------------------------------------------------

sim_results_output$setting_table_ssd <- setting_table_ssd
sim_results_output$setting_table_any <- setting_table_any

save(sim_results_output, file = paste0(sim_out_path,"sim_results_output.RData"))


############################################
#### Create Readme of analysis overview ####
############################################

# add readme noting build info
heading <- paste0("Simulation for ",cond_name," using any visits")
sim_date <- paste0("Simulation Date: ",Sys.time())
trial_size <- paste0("Simulation Trials: ",n_trials)
population_size <- paste0("Total Population: ",n_trials)
change_point_used <- paste0("Change-point used: ",delay_params$cp," days")
optimal_model <- paste0("Selected Model: ",selected_model$model)
tmp <- sim_data_all$miss_bins_visits %>% summarise(n_miss =sum(num_miss)) %>% .$n_miss
total_miss <- paste0("Total Number of misses: ",tmp)
tmp1 <- names(other_model_res)
tmp2 <- map(other_model_res,~.$total_miss)
total_miss_other <- c("Total miss other models: \n",paste0("  ",tmp1," - ",tmp2," \n"))

# collect results
tmp <- map(other_model_res,~.$model_res$sim_res_stats)
tmp[[selected_model$model]] <- optimal_model_res$sim_res_stats

pct_miss_summary <- c("Percent Miss: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[1])," \n"))
mean_dur_summary <- c("Mean Duration: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[3])," \n"))
median_dur_summary <- c("Median Duration: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[4])," \n"))
mean_n_miss_summary <- c("Mean Number of Misses (per patient): \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[5])," \n"))
median_n_miss_summary <- c("Median Number of Misses (per patient): \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[6])," \n"))


fileConn<-file(paste0(sim_out_path,"readme.txt"))
writeLines(c(heading,
             sim_date,
             trial_size,
             population_size,
             change_point_used,
             optimal_model,
             total_miss,
             total_miss_other,
             pct_miss_summary,
             mean_dur_summary,
             median_dur_summary,
             mean_n_miss_summary,
             median_n_miss_summary), fileConn)
close(fileConn)





