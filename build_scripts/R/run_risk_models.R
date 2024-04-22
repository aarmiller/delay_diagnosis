
library(tidyverse)
library(lubridate)

args = commandArgs(trailingOnly=TRUE)
# name of condition
cond_name <- args[1]
# cond_name <- "dengue"

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

delay_base_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")
# delay_base_path <- "/Volumes/AML/small_dbs/tb/truven/enroll_restrict_365/delay_results/"

out_path <- paste0(delay_base_path,"risk_models/")

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

## Setup model parameters ------------------------------------------------------
num_cores <- 40

## Setup age categories
age_cats <- c(-1,17,34,44,54,64,130)

## Main setting labels
setting_labels <- tribble(~outpatient,~ed,~obs_stay,~inpatient,~setting_label,
                          1,0,0,0,"Out only",
                          1,1,0,0,"Out and ED",
                          1,0,1,0,"Out and Obs",
                          1,0,0,1,"Out and Inpatient",
                          1,1,1,0,"Out and ED and Obs",
                          1,1,0,1,"Out and ED and Inpatient",
                          1,0,1,1,"Out and Obs and Inpatient",
                          1,1,1,1,"All Four Settings",
                          0,1,0,0,"ED only",
                          0,1,1,0,"ED and Obs",
                          0,1,0,1,"ED and Inpatient",
                          0,1,1,1,"ED and Obs and Inpatient",
                          0,0,1,0,"Obs only",
                          0,0,1,1,"Obs and Inpatient",
                          0,0,0,1,"Inpatient only") %>% 
  mutate(setting_label = fct_relevel(setting_label,
                                     "Out only",
                                     "ED only",
                                     "Obs only",
                                     "Inpatient only",
                                     "Out and ED",
                                     "Out and Obs",
                                     "Out and Inpatient",
                                     "ED and Obs",
                                     "ED and Inpatient",
                                     "Obs and Inpatient",
                                     "Out and ED and Obs",
                                     "Out and ED and Inpatient",
                                     "Out and Obs and Inpatient",
                                     "ED and Obs and Inpatient",
                                     "All Four Settings")) 


## Load Delay Data -------------------------------------------------------------

load(paste0(delay_base_path,"demo_data.RData"))
# demo1
# demo2
# rural_visits

load(paste0(delay_base_path,"all_dx_visits.RData"))
# visit_counts
# all_dx_visits

load(paste0(delay_base_path,"delay_tm.RData"))
# tm

load(paste0(delay_base_path,"ssd_visit/sim_res.RData"))
# sim_res
# sim_res_sim_obs

n_trials <- max(sim_res$trial)



############################
#### Prepare Visit Info ####
############################


index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct(patient_id,outpatient,ed,obs_stay,inpatient,svcdate) %>% 
  mutate(dow=weekdays(as_date(svcdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow) %>%
  inner_join(setting_labels,by = join_by(outpatient, ed, obs_stay, inpatient))


obs_locations <- tm %>% 
  left_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct(obs,patient_id,outpatient,ed,obs_stay,inpatient,svcdate) %>%
  inner_join(setting_labels,by = join_by(outpatient, ed, obs_stay, inpatient))


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

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))


### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res <- sim_res %>% 
  select(obs,trial) %>% 
  inner_join(obs_locations,by = join_by(obs)) %>% 
  inner_join(sim_res_sim_obs %>% 
               inner_join(select(index_dx_dates,patient_id,index_date), by = "patient_id") %>% 
               mutate(vis_date = index_date+days_since_index) %>% 
               mutate(dow=weekdays(as_date(vis_date))) %>% 
               select(obs,dow) %>% 
               mutate(weekend = dow %in% c("Saturday","Sunday")), by = "obs")



###########################
#### Regression Models ####
###########################

#### Missed opportunities All --------------------------------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id")
  
  
  fit <- glm(miss~setting_label+obs_stay + age_cat + female + rural + source + weekend, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_miss_res(10)

cluster <- parallel::makeCluster(num_cores)
# cluster <- parallel::makeCluster(6)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("sim_res","get_miss_res","index_locations","reg_demo"),
                        envir=environment())


miss_opp_res <- parallel::parLapply(cl = cluster,
                                    1:max(sim_res$trial),
                                    function(x){get_miss_res(x)})


# parallel::stopCluster(cluster)

miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate),na.rm = TRUE),
            low = quantile(exp(estimate),probs = c(0.025),na.rm = TRUE),
            high = quantile(exp(estimate),probs = c(0.975),na.rm = TRUE))

parallel::stopCluster(cluster)
rm(cluster)
gc()


## Using 2 (inpatient vs outpatient) setting labels ----------------------------

get_miss_res2 <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id")
  
  
  fit <- glm(miss~inpatient+obs_stay + age_cat + female + rural + source + weekend, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

cluster <- parallel::makeCluster(num_cores)
# cluster <- parallel::makeCluster(6)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("sim_res","get_miss_res2","index_locations","reg_demo"),
                        envir=environment())


miss_opp_res2 <- parallel::parLapply(cl = cluster,
                                     1:max(sim_res$trial),
                                     function(x){get_miss_res2(x)})


# parallel::stopCluster(cluster)

miss_opp_res2 <- bind_rows(miss_opp_res2) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate),na.rm = TRUE),
            low = quantile(exp(estimate),probs = c(0.025),na.rm = TRUE),
            high = quantile(exp(estimate),probs = c(0.975),na.rm = TRUE))

parallel::stopCluster(cluster)
rm(cluster)
gc()


#### Missed opportunities Medicaid ---------------------------------------------

get_miss_res_med <- function(trial_val){
  
  tmp1 <- sim_res %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(miss~setting_label+obs_stay + age_cat + female + rural + race + weekend, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_miss_res_med(10)

cluster <- parallel::makeCluster(num_cores)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("get_miss_res_med","sim_res","index_locations","reg_demo"),
                        envir=environment())


miss_opp_res_med <- parallel::parLapply(cl = cluster,
                                        1:max(sim_res$trial),
                                        function(x){get_miss_res_med(x)})


parallel::stopCluster(cluster)
rm(cluster)
gc()

miss_opp_res_med <- bind_rows(miss_opp_res_med) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate),na.rm = TRUE),
            low = quantile(exp(estimate),probs = c(0.025),na.rm = TRUE),
            high = quantile(exp(estimate),probs = c(0.975),na.rm = TRUE))


## Using 2 (inpatient vs outpatient) setting labels ----------------------------

get_miss_res_med2 <- function(trial_val){
  
  tmp1 <- sim_res %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(miss~inpatient+obs_stay + age_cat + female + rural + race + weekend, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_miss_res_med(10)

cluster <- parallel::makeCluster(num_cores)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("get_miss_res_med2","sim_res","index_locations","reg_demo"),
                        envir=environment())


miss_opp_res_med2 <- parallel::parLapply(cl = cluster,
                                         1:max(sim_res$trial),
                                         function(x){get_miss_res_med2(x)})


parallel::stopCluster(cluster)
rm(cluster)
gc()

miss_opp_res_med2 <- bind_rows(miss_opp_res_med2) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate),na.rm = TRUE),
            low = quantile(exp(estimate),probs = c(0.025),na.rm = TRUE),
            high = quantile(exp(estimate),probs = c(0.975),na.rm = TRUE))



#### Delay duration All --------------------------------------------------------

# compute duration by simulation (note this will be used in the next step)
sim_res_dur <- sim_res %>% 
  select(obs,trial) %>% 
  inner_join(sim_res_sim_obs,by = "obs") %>% 
  group_by(trial,patient_id) %>% 
  summarise(duration = -min(days_since_index)) %>% 
  ungroup()

rm(sim_res)
gc()

get_dur_res <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L))
  
  
  fit <- glm(duration~age_cat + female + rural + source, family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_dur_res(10)

cluster <- parallel::makeCluster(num_cores)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("sim_res_dur","get_dur_res","reg_demo"),
                        envir=environment())


miss_dur_res <- parallel::parLapply(cl = cluster,
                                    1:n_trials,
                                    function(x){get_dur_res(x)})


# parallel::stopCluster(cluster)

miss_dur_res <- bind_rows(miss_dur_res) %>% 
  group_by(term) %>% 
  summarise(est = mean((estimate)),
            low = quantile((estimate),probs = c(0.025)),
            high = quantile((estimate),probs = c(0.975)))


#### Delay duration Medicaid ---------------------------------------------------



get_dur_res_med <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(duration~age_cat + female + rural + race, family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_dur_res_med(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("get_dur_res_med"),
                        envir=environment())


miss_dur_res_med <- parallel::parLapply(cl = cluster,
                                        1:n_trials,
                                        function(x){get_dur_res_med(x)})

# parallel::stopCluster(cluster)

miss_dur_res_med <- bind_rows(miss_dur_res_med) %>% 
  group_by(term) %>% 
  summarise(est = mean((estimate)),
            low = quantile((estimate),probs = c(0.025)),
            high = quantile((estimate),probs = c(0.975)))




#### Delay Patient All ---------------------------------------------------------


get_delay_pat_res <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    mutate(miss = duration>0)
  
  
  fit <- glm(miss~age_cat + female + rural + source, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_delay_pat_res(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("get_delay_pat_res"),
                        envir=environment())


miss_delay_pat_res <- parallel::parLapply(cl = cluster,
                                          1:n_trials,
                                          function(x){get_delay_pat_res(x)})


# parallel::stopCluster(cluster)

miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate)),
            low = quantile(exp(estimate),probs = c(0.025)),
            high = quantile(exp(estimate),probs = c(0.975)))


#### Delay Patient Medicaid ----------------------------------------------------


get_delay_pat_res_med <- function(trial_val){
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    mutate(miss = duration>0) %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(miss~age_cat + female + rural + race, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_delay_pat_res_med(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("get_delay_pat_res_med"),
                        envir=environment())


miss_delay_pat_res_med <- parallel::parLapply(cl = cluster,
                                              1:n_trials,
                                              function(x){get_delay_pat_res_med(x)})


parallel::stopCluster(cluster)
rm(cluster)
gc()

miss_delay_pat_res_med <- bind_rows(miss_delay_pat_res_med) %>% 
  group_by(term) %>% 
  summarise(est = mean(exp(estimate)),
            low = quantile(exp(estimate),probs = c(0.025)),
            high = quantile(exp(estimate),probs = c(0.975)))





### Alternative missed visits --------------------------------------------------
# 
# obs_locations2 <- tm %>% 
#   left_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
#   filter(!is.na(obs)) %>% 
#   filter(setting_type !=4) %>%
#   distinct(patient_id,obs,setting_type,admdate) %>% 
#   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
#   mutate(dow=weekdays(as_date(admdate))) %>% 
#   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
#   select(patient_id,obs,setting,dow,weekend) 
# 
# index_locations2 <- tm %>% 
#   filter(days_since_index==0) %>% 
#   filter(setting_type !=4) %>%
#   distinct(patient_id,setting_type,admdate) %>% 
#   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
#   mutate(dow=weekdays(as_date(admdate))) %>% 
#   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
#   select(patient_id,setting,weekend,dow) 
# 
# sim_res <- sim_res %>% 
#   select(obs,trial) %>% 
#   inner_join(obs_locations2)
# 
# get_miss_res2 <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- sim_res %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% 
#                           mutate(miss=TRUE),
#                         index_locations2 %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo, by = "patient_id")
#   
#   
#   fit <- glm(miss~setting + age_cat + female + rural + source + weekend, family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# get_miss_res2(10)
# 
# cluster <- parallel::makeCluster(num_cores)
# # cluster <- parallel::makeCluster(6)
# 
# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("sim_res","get_miss_res2","index_locations2","reg_demo"),
#                         envir=environment())
# 
# 
# miss_opp_res2 <- parallel::parLapply(cl = cluster,
#                                     1:max(sim_res$trial),
#                                     function(x){get_miss_res2(x)})
# 
# 
# parallel::stopCluster(cluster)
# 
# miss_opp_res2 <- bind_rows(miss_opp_res2) %>% 
#   group_by(term) %>% 
#   summarise(est = mean(exp(estimate)),
#             low = quantile(exp(estimate),probs = c(0.025)),
#             high = quantile(exp(estimate),probs = c(0.975)))
# 
# 
# miss_opp_res2
# miss_opp_res

########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res,
                             miss_opp_res2 = miss_opp_res2,
                             miss_opp_res_med = miss_opp_res_med,
                             miss_opp_res_med2 = miss_opp_res_med2,
                             miss_dur_res = miss_dur_res,
                             miss_dur_res_med = miss_dur_res_med,
                             miss_delay_pat_res = miss_delay_pat_res,
                             miss_delay_pat_res_med = miss_delay_pat_res_med)


save(ssd_miss_risk_models,file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = "github/delay_diagnosis/build_scripts/R/report_scripts/risk_model_report.Rmd",
                  params = list(cond = cond_name),
                  output_dir = out_path)