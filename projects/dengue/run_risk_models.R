
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "dengue"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

num_cores <- 40
cp_selected <- delay_params$cp

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
tm <- tm %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

load(paste0(sim_out_path,"/sim_res_ssd.RData"))
# sim_res_ssd


n_trials <- nrow(distinct(sim_res_ssd))

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


## Prepare antibiotic indicators -----------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

rx_visits <- con %>% tbl("all_rx_visits") %>% 
  filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
  collect()

rx_visits <- rx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))

gc()

abx_codes <- read_csv(paste0(delay_params$out_path, "antibiotics.csv"))
abx_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% abx_codes$ndcnum) 

reg_demo <- reg_demo %>% 
  left_join(abx_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(abx_window=1L)) %>% 
  mutate_at(vars(abx_window),~replace_na(.,0L))

save(reg_demo, file=paste0(out_path,"reg_data.RData"))
## Prepare ID consult indicators -----------------------------------------------

# Aaron wants to know among those who received testing for dengue 
# how close was testing date to date of consult
# proc_codes <- c("86790", "87449", "87798")
# stdproc_codes <- c("285", "448")
# procs <- con %>%
#   tbl("all_proc_visits") %>%
#   filter(proc %in% proc_codes) %>%
#   filter(between(days_since_index,-14,14)) %>%
#   collect() %>% inner_join(tbl(con, "index_dx_dates") %>% collect() %>%
#                              select(patient_id,index_date )) %>%
#   mutate(test_date = index_date+days_since_index)
# 
# index_case_date <- tbl(con, "index_dx_dates") %>% collect() %>%
#   inner_join(index_cases)
# 
# plan <- collect_plan(con)
# DBI::dbListTables(con)
# 
# temp.out <- plan %>% mutate(setting = map(year, ~c("facility", "outpatient"))) %>% unnest(setting) %>%
#   filter(!(year == "01" & setting == "facility")) %>% 
#   gether_core_data(vars = c("patient_id", "stdplac", "svcdate",
#                             "stdprov"),
#                    db_con = con) %>% dplyr::mutate(core_data = purrr::map(core_data,
#                                                                              ~dplyr::distinct(.)))
# temp.out <- temp.out %>% dplyr::mutate(source_type = ifelse(source ==
#                                                               "ccae", 1L, ifelse(source == "mdcr", 0L, 2L))) %>% dplyr::select("year",
#                                                                                                                                "source_type", "core_data") %>% mutate(core_data = map(core_data,
#                                                                                                                                                                                       ~mutate(., stdprov = as.character(stdprov)))) %>%
#   tidyr::unnest(cols = c(core_data)) %>% dplyr::mutate(disdate = svcdate,
#                                                        admdate = svcdate,) %>% dplyr::select(-svcdate) %>%
#   distinct()
# 
# temp.out <- temp.out %>% inner_join(index_case_date) %>% distinct() %>%
#   filter(stdprov %in% stdproc_codes) %>%
#   filter(admdate<=index_date+14 & (admdate)>=(index_date-14)) %>% # ID consult within +/-14 days of index
#   distinct(index_date,admdate, patient_id)
# 
# # among those that recieved testing for dengue, could have multiple ID consults
# inner_join(procs, temp.out) %>%
#   mutate(time_between_test_ID = test_date-admdate) %>%
#   count(time_between_test_ID) %>% print(n = Inf)

# This is the old function to find consult ids before Aaron added the stdprov_visits table to the database
# consult_ids <- function (db_con, index_cases, stdproc_codes) {
#   plan <- collect_plan(db_con)
#   temp.out <- plan %>% mutate(setting = map(year, ~c("facility", "outpatient"))) %>% unnest(setting) %>%
#     filter(!(year == "01" & setting == "facility")) %>% 
#     gether_core_data(vars = c("patient_id", "stdplac", "svcdate",
#                               "stdprov"),
#                      db_con = db_con) %>% dplyr::mutate(core_data = purrr::map(core_data,
#                                                                                ~dplyr::distinct(.)))
#   temp.out <- temp.out %>% dplyr::mutate(source_type = ifelse(source ==
#                                                                 "ccae", 1L, ifelse(source == "mdcr", 0L, 2L))) %>% dplyr::select("year",
#                                                                                                                                  "source_type", "core_data") %>% mutate(core_data = map(core_data,
#                                                                                                                                                                                         ~mutate(., stdprov = as.character(stdprov)))) %>%
#     tidyr::unnest(cols = c(core_data)) %>% dplyr::mutate(disdate = svcdate,
#                                                          admdate = svcdate) %>% dplyr::select(-svcdate) %>%
#     distinct()
#   
#   temp.out <- temp.out %>% inner_join(index_cases, by = "patient_id") %>% distinct() %>% 
#     filter(stdprov %in% stdproc_codes) %>% 
#     filter(admdate<=index_date & (admdate)>=(index_date-(delay_params$upper_bound))) %>% 
#     distinct(patient_id, admdate)
#   return(temp.out)
# }
# 
# 
# ID_consult_vis_old <- consult_ids(con, 
#                               index_cases %>% 
#                                 mutate(index_date = index_date + shift) %>% 
#                                 select(patient_id, index_date),
#                               c("285", "448")) 


ID_consult_vis <- con %>% tbl("stdprov_visits") %>% 
  filter(patient_id %in% local(index_cases$patient_id)) %>% 
  filter(stdprov %in%  c("285", "448")) %>% 
  distinct() %>% 
  collect() %>% 
  inner_join(index_cases %>% 
             mutate(index_date = index_date + shift) %>% 
             select(patient_id, index_date), by = "patient_id") %>% 
  rename(admdate = svcdate) %>% 
  filter(admdate<=index_date & (admdate)>=(index_date-(delay_params$upper_bound))) %>% 
  distinct(patient_id, admdate)

## Add abx and ID consult to time map
abx_rx_vis <- abx_rx_vis %>% distinct(patient_id, admdate = date)

tm <- tm %>% left_join(abx_rx_vis %>% mutate(abx = 1L), by = c("patient_id", "admdate")) %>% 
  left_join(ID_consult_vis %>% mutate(ID_consult = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(abx, ID_consult),~replace_na(.,0L)) 

############################
#### Prepare Visit Info ####
############################

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,setting_type,admdate, abx, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  select(patient_id,setting,weekend,dow, abx, ID_consult) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(tribble(~outpatient,~ed,~inpatient,~setting_label,
                     T,F,F,"Out only",
                     T,T,F,"Out and ED",
                     T,T,T,"All three",
                     T,F,T,"Out and inpatient",
                     F,T,F,"ED only",
                     F,T,T,"ED and inpatient",
                     F,F,T,"Inpatient only",
                     F,F,F,"none")) %>% 
  mutate(setting_label = fct_relevel(setting_label,
                                     "Out only",
                                     "ED only",
                                     "Inpatient only",
                                     "Out and ED",
                                     "Out and inpatient",
                                     "ED and inpatient",
                                     "All three",
                                     "none")) 

obs_locations <- tm %>% 
  left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,obs,setting_type, abx, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  select(patient_id,obs,setting, abx, ID_consult) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(tribble(~outpatient,~ed,~inpatient,~setting_label,
                     T,F,F,"Out only",
                     T,T,F,"Out and ED",
                     T,T,T,"All three",
                     T,F,T,"Out and inpatient",
                     F,T,F,"ED only",
                     F,T,T,"ED and inpatient",
                     F,F,T,"Inpatient only",
                     F,F,F,"none")) %>% 
  mutate(setting_label = fct_relevel(setting_label,
                                     "Out only",
                                     "ED only",
                                     "Inpatient only",
                                     "Out and ED",
                                     "Out and inpatient",
                                     "ED and inpatient",
                                     "All three",
                                     "none")) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd <- sim_res_ssd %>% 
  select(res) %>% 
  mutate(trial = row_number()) %>% 
  unnest(res) %>% 
  inner_join(obs_locations) %>% 
  inner_join(sim_obs %>% 
               inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
                            select(patient_id,index_date), by = "patient_id") %>% 
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
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id")
  
  
  fit <- glm(miss~setting_label + obs_stay + age_cat + female + rural + source + weekend + abx + ID_consult, 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_miss_res(10)

# cluster <- parallel::makeCluster(num_cores)
# # cluster <- parallel::makeCluster(6)
# 
# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("sim_res_ssd","get_miss_res","index_locations","reg_demo"),
#                         envir=environment())
# 
# 
# miss_opp_res <- parallel::parLapply(cl = cluster,
#                                     1:max(sim_res_ssd$trial),
#                                     function(x){get_miss_res(x)})
# parallel::stopCluster(cluster)

miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                            function(x){get_miss_res(x)}, 
                                            mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# rm(cluster)
gc()


#### Missed opportunities Medicaid ---------------------------------------------

get_miss_res_med <- function(trial_val){
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(miss~setting_label+obs_stay + age_cat + female + rural + race + weekend + abx + ID_consult,
             family = "binomial", data=reg_data)
  
  # fit.penalized <- logistf(miss~setting_label+obs_stay + age_cat + female + rural + race + weekend + abx ,
  #                          family = binomial, data=reg_data,
  #                          control = logistf.control( maxit = 1000),
  #                          plcontrol = logistpl.control( maxit = 1000)) 
  # 
  # tibble(term = names(coef(fit.penalized)),
  #        estimate = coef(fit.penalized))
  
  broom::tidy(fit)
  
  
}

# get_miss_res_med(10)

# cluster <- parallel::makeCluster(num_cores)
# 
# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("get_miss_res_med","sim_res_ssd","index_locations","reg_demo"),
#                         envir=environment())
# 
# 
# miss_opp_res_med <- parallel::parLapply(cl = cluster,
#                                         1:max(sim_res_ssd$trial),
#                                         function(x){get_miss_res_med(x)})
# 
# 
# parallel::stopCluster(cluster)


miss_opp_res_med <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                            function(x){get_miss_res_med(x)}, 
                                            mc.cores = num_cores)

miss_opp_res_med <- bind_rows(miss_opp_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


# rm(cluster)
gc()


#### Delay duration All --------------------------------------------------------

# compute duration by simulation (note this will be used in the next step)
sim_res_dur <- sim_res_ssd %>% 
  select(obs,trial) %>% 
  inner_join(sim_obs,by = "obs") %>% 
  group_by(trial,patient_id) %>% 
  summarise(duration = -min(days_since_index)) %>% 
  ungroup()

rm(sim_res_ssd)
gc()

get_dur_res <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L))
  
  
  fit <- glm(duration~age_cat + female + rural + source + abx_window, family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_dur_res(10)

# cluster <- parallel::makeCluster(num_cores)
# 
# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("sim_res_dur","get_dur_res","reg_demo"),
#                         envir=environment())
# 
# 
# miss_dur_res <- parallel::parLapply(cl = cluster,
#                                     1:n_trials,
#                                     function(x){get_dur_res(x)})
# 
# parallel::stopCluster(cluster)

miss_dur_res <- parallel::mclapply(1:max(sim_res_dur$trial),
                                       function(x){get_dur_res(x)}, 
                                       mc.cores = num_cores)


miss_dur_res <- bind_rows(miss_dur_res) %>% 
  group_by(term) %>% 
  summarise(est = median(estimate, na.rm = T),
            low = quantile(estimate, probs = c(0.025), na.rm = T),
            high = quantile(estimate, probs = c(0.975), na.rm = T))


gc()

#### Delay duration Medicaid ---------------------------------------------------

get_dur_res_med <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(duration~age_cat + female + rural + race + abx_window, family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_dur_res_med(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("get_dur_res_med"),
#                         envir=environment())
# 
# 
# miss_dur_res_med <- parallel::parLapply(cl = cluster,
#                                         1:n_trials,
#                                         function(x){get_dur_res_med(x)})
# 
# parallel::stopCluster(cluster)

miss_dur_res_med <- parallel::mclapply(1:max(sim_res_dur$trial),
                                   function(x){get_dur_res_med(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_med <- bind_rows(miss_dur_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(estimate, na.rm = T),
            low = quantile(estimate ,probs = c(0.025), na.rm = T),
            high = quantile(estimate ,probs = c(0.975), na.rm = T))

gc()

#### Delay Patient All ---------------------------------------------------------


get_delay_pat_res <- function(trial_val){
  
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    mutate(miss = duration>0)
  
  
  fit <- glm(miss~age_cat + female + rural + source + abx_window, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_delay_pat_res(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("get_delay_pat_res"),
#                         envir=environment())
# 
# 
# miss_delay_pat_res <- parallel::parLapply(cl = cluster,
#                                           1:n_trials,
#                                           function(x){get_delay_pat_res(x)})

# parallel::stopCluster(cluster)

miss_delay_pat_res <- parallel::mclapply(1:max(sim_res_dur$trial),
                                       function(x){get_delay_pat_res(x)}, 
                                       mc.cores = num_cores)

miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

gc()

#### Delay Patient Medicaid ----------------------------------------------------


get_delay_pat_res_med <- function(trial_val){
  tmp1 <- sim_res_dur %>% filter(trial==trial_val)
  
  reg_data <- reg_demo %>% 
    left_join(tmp1) %>% 
    mutate(duration = replace_na(duration,0L)) %>% 
    mutate(miss = duration>0) %>% 
    filter(source == "medicaid")
  
  
  fit <- glm(miss~age_cat + female + rural + race + abx_window, family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}

# get_delay_pat_res_med(10)

# cluster <- parallel::makeCluster(num_cores)

# parallel::clusterCall(cluster, function() library(tidyverse))
# parallel::clusterExport(cluster,c("get_delay_pat_res_med"),
#                         envir=environment())
# 
# 
# miss_delay_pat_res_med <- parallel::parLapply(cl = cluster,
#                                               1:n_trials,
#                                               function(x){get_delay_pat_res_med(x)})
# 
# 
# parallel::stopCluster(cluster)

miss_delay_pat_res_med <- parallel::mclapply(1:max(sim_res_dur$trial),
                                         function(x){get_delay_pat_res_med(x)}, 
                                         mc.cores = num_cores)

miss_delay_pat_res_med <- bind_rows(miss_delay_pat_res_med) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


gc()

### Alternative missed visits --------------------------------------------------
# 
# obs_locations2 <- tm %>% 
#   left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
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
# sim_res_ssd <- sim_res_ssd %>% 
#   select(obs,trial) %>% 
#   inner_join(obs_locations2)
# 
# get_miss_res2 <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- sim_res_ssd %>% 
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
# parallel::clusterExport(cluster,c("sim_res_ssd","get_miss_res2","index_locations2","reg_demo"),
#                         envir=environment())
# 
# 
# miss_opp_res2 <- parallel::parLapply(cl = cluster,
#                                     1:max(sim_res_ssd$trial),
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
                             miss_opp_res_med = miss_opp_res_med,
                             miss_dur_res = miss_dur_res,
                             miss_dur_res_med = miss_dur_res_med,
                             miss_delay_pat_res = miss_delay_pat_res,
                             miss_delay_pat_res_med = miss_delay_pat_res_med)


save(ssd_miss_risk_models,file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/dengue/risk_model_report.Rmd",
                  params = list(cond = "Dengue"),
                  output_dir = out_path)