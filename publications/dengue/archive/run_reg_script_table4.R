
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

out_path <-paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/")

### Load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

### Connect to db
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))


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

# problem patient_ids that do not have index visits in tm (i.e. days_since_index ==0)
problem_patient_ids <- tm %>% filter(days_since_index == 0) %>% distinct(patient_id) %>%
  anti_join(tm %>% distinct(patient_id),.)

# Get setting type for index dates
setting_problem_patient_id <- tbl(con, "tm_full") %>% filter(patient_id %in% local(problem_patient_ids$patient_id)) %>% 
  collect() %>% inner_join(index_cases %>% select(patient_id, svcdate = index_date) %>% 
                             inner_join(problem_patient_ids)) %>% 
  select(patient_id, admdate = svcdate, disdate, outpatient:inpatient) %>% 
  pivot_longer(outpatient:inpatient, names_to = "setting_type", values_to = "ind") %>% 
  filter(ind !=0) %>% select(-ind) %>% 
  mutate(setting_type = as.integer(factor(setting_type, labels= 1:5,
                                          levels = c("outpatient",
                                                     "ed",
                                                     "obs_stay",
                                                     "rx",
                                                     "inpatient"))))

# Get stdplac type for index dates
std_plac_problem_patient_id <- tbl(con, "stdplac_visits") %>% filter(patient_id %in% local(problem_patient_ids$patient_id)) %>% 
  collect() %>% inner_join(index_cases %>% select(patient_id, svcdate = index_date) %>% 
                             inner_join(problem_patient_ids)) %>% 
  select(patient_id, admdate = svcdate, stdplac)

# construct tm for index dates for problem vis
tm_problem_patient_id <- setting_problem_patient_id %>% inner_join(std_plac_problem_patient_id) %>% 
  inner_join(index_cases %>% select(patient_id, index_date) %>% 
               inner_join(problem_patient_ids)) %>% 
  mutate(days_since_index = admdate - index_date) %>% 
  select(names(tm))

tm <- tm %>% bind_rows(tm_problem_patient_id) %>% 
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

opioid_codes <- codeBuildr::load_rx_codes("opioids") %>% unlist(use.names = F) %>% unique()

opioid_rx_vis <- rx_visits %>% 
  filter(ndcnum %in% opioid_codes) 

reg_demo <- reg_demo %>% 
  left_join(opioid_rx_vis %>% 
              filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
              distinct(patient_id) %>% 
              mutate(opioid_window=1L)) %>% 
  mutate_at(vars(opioid_window),~replace_na(.,0L))

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
opioid_rx_vis <- opioid_rx_vis %>% distinct(patient_id, admdate = date)


tm <- tm %>% left_join(abx_rx_vis %>% mutate(abx = 1L), by = c("patient_id", "admdate")) %>% 
  left_join(opioid_rx_vis %>% mutate(opioid = 1L), by = c("patient_id", "admdate")) %>% 
  left_join(ID_consult_vis %>% mutate(ID_consult = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(abx, opioid, ID_consult),~replace_na(.,0L)) 

############################
#### Prepare Visit Info ####
############################

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,setting_type,admdate, abx, opioid, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  select(patient_id,setting,weekend,dow, year, month, abx, opioid, ID_consult) %>%
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
  distinct(patient_id,obs,setting_type, abx, opioid, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  select(patient_id,obs,setting, abx, opioid, ID_consult) %>%
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


setting_labels <- expand.grid(outpatient = c("Outpatient", NA),
                              inpatient = c("Inpatient", NA),
                              ed = c("ED", NA),
                              obs_stay = c("Observational Stay", NA)) %>% 
  unite(., col = setting_label, sep = " + ", na.rm = T, remove =F) %>% 
  mutate(across(outpatient:obs_stay, ~ifelse(is.na(.), F, T))) %>% 
  mutate(setting_label = ifelse(setting_label == "", "Not any",
                                ifelse(setting_label == "Inpatient", "Inpatient Only",
                                       ifelse(setting_label == "Outpatient", "Outpatient Only",
                                              ifelse(setting_label == "ED", "ED Only",
                                                     ifelse(setting_label == "Observational Stay", "Observational Stay Only",
                                                            setting_label))))))
setting_labels <- setting_labels[order(rowSums(setting_labels[,2:5])), ]
setting_labels <- setting_labels %>% mutate(setting_label = factor(setting_label, 
                                                                   levels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                                                              setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                                                              setting_labels$setting_label[which(setting_labels$setting_label=="Not any")]),
                                                                   labels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                                                              setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                                                              setting_labels$setting_label[which(setting_labels$setting_label=="Not any")])))

index_locations_2 <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,setting_type,admdate, abx, opioid, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  select(patient_id,setting,weekend,dow, year, month, abx, opioid, ID_consult) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(setting_labels) 

obs_locations_2 <- tm %>% 
  left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,obs,setting_type, abx, opioid, ID_consult) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  select(patient_id,obs,setting, abx, opioid, ID_consult) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(setting_labels) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd_2 <- sim_res_ssd %>% 
  select(res) %>% 
  mutate(trial = row_number()) %>% 
  unnest(res) %>% 
  inner_join(obs_locations_2) %>% 
  inner_join(sim_obs %>% 
               inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
                            select(patient_id,index_date), by = "patient_id") %>% 
               mutate(vis_date = index_date+days_since_index) %>% 
               mutate(year = year(as_date(vis_date)),
                      month = month(as_date(vis_date))) %>% 
               mutate(dow=weekdays(as_date(vis_date))) %>% 
               select(obs, dow, year, month) %>% 
               mutate(weekend = dow %in% c("Saturday","Sunday")), by = "obs")

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
               mutate(year = year(as_date(vis_date)),
                      month = month(as_date(vis_date))) %>% 
               mutate(dow=weekdays(as_date(vis_date))) %>% 
               select(obs, dow, year, month) %>% 
               mutate(weekend = dow %in% c("Saturday","Sunday")), by = "obs")

out_path <-paste0("/Shared/Statepi_Diagnosis/projects/", proj_name, "/risk_models/")

###########################
#### Regression Models ####
###########################


#### Missed opportunities All (original setting labels)-------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>%
                          mutate(miss=TRUE),
                        index_locations %>%
                          mutate(miss=FALSE)) %>%
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                          levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
    
  fit <- glm(miss~setting_label + obs_stay + age_cat + female + source + weekend + abx + ID_consult
             # + year + month
             , 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)

gc()

original_setting_labels <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# table4 <- table4 %>% 
#   mutate(nice = paste0(trimws(format(round(est, 2), nsmall = 2)),
#                        " (", 
#                        trimws(format(round(low, 2), nsmall = 2)),
#                        "-",
#                        trimws(format(round(high, 2), nsmall = 2)),
#                        ")"))
# rm(cluster)
# write_csv(table4, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table4.csv"))
# 

# save(original_setting_labels, file = paste0(out_path, "reg_data/original_setting_labels.RData"))
# save(original_setting_labels, file = paste0(out_path, "reg_data/original_setting_labels_g_2001.RData"))
save(original_setting_labels, file = paste0(out_path, "reg_data/original_setting_labels_no_year_month.RData"))


#### Missed opportunities All (original setting labels - Firth regression)-------------------------
library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>%
                          mutate(miss=TRUE),
                        index_locations %>%
                          mutate(miss=FALSE)) %>%
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  # fit <- glm(miss~setting_label + obs_stay + age_cat + female + source + weekend + abx + ID_consult + year + month, 
  #            family = "binomial", data=reg_data)
  
  fit <- logistf(miss~setting_label + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
                 # + year + month
                 ,
                 family = "binomial", data=reg_data,
                 control = logistf.control( maxit = 1000),
                 plcontrol = logistpl.control( maxit = 1000))

  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  # broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

original_setting_labels_firth <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# save(original_setting_labels_firth, file = paste0(out_path, "reg_data/original_setting_labels_firth.RData"))
# save(original_setting_labels_firth, file = paste0(out_path, "reg_data/original_setting_labels_firth_g_2001.RData"))
save(original_setting_labels_firth, file = paste0(out_path, "reg_data/original_setting_labels_firth_no_year_month.RData"))


#### Missed opportunities All (original setting labels - remove sep)-------------------------
library(detectseparation)

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  
  test <- glm(miss~setting_label + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
              # + year + month
              ,
              family = "binomial", data=reg_data, method = "detect_separation")
  
  if(!test$outcome){
    fit <- glm(miss~setting_label + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
               # + year + month
               ,
               family = "binomial", data=reg_data)
    
    out <- tibble(term = names(coef(fit)),
                  estimate = coef(fit))
  } else {
    out <- tibble(term = "quasi_sep",
                  estimate = NA)
  }
  return(out)
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

original_setting_labels_remove <- bind_rows(miss_opp_res) %>% 
  filter(term != "quasi_sep") %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

count_quasi_sep <- bind_rows(miss_opp_res) %>% 
  filter(term == "quasi_sep") %>% 
  count() %>% .$n

original_setting_labels_remove <- bind_rows(original_setting_labels_remove,
                                            tibble(term = "sim_trials_removed",
                                                   est = count_quasi_sep,
                                                   low = NA,
                                                   high = NA))
  
# save(original_setting_labels_remove, file = paste0(out_path, "reg_data/original_setting_labels_remove.RData"))
# save(original_setting_labels_remove, file = paste0(out_path, "reg_data/original_setting_labels_remove_g_2001.RData"))
save(original_setting_labels_remove, file = paste0(out_path, "reg_data/original_setting_labels_remove_no_year_month.RData"))


################################################################################

#### Missed opportunities All (new setting labels)-------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd_2 %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations_2 %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  fit <- glm(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult
             # + year + month
             , 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

new_setting_labels <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


# save(new_setting_labels, file = paste0(out_path, "reg_data/new_setting_labels.RData"))
# save(new_setting_labels, file = paste0(out_path, "reg_data/new_setting_labels_g_2001.RData"))
save(new_setting_labels, file = paste0(out_path, "reg_data/new_setting_labels_no_year_month.RData"))



#### Missed opportunities All (new setting labels - Firth regression)-------------------------
library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd_2 %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations_2 %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
   
  # fit <- glm(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult + year + month, 
  #            family = "binomial", data=reg_data)
  
  fit <- logistf(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult 
                 # + year + month
                 ,
                 family = "binomial", data=reg_data,
                 control = logistf.control( maxit = 1000),
                 plcontrol = logistpl.control( maxit = 1000))
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  # broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

new_setting_labels_firth <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# save(new_setting_labels_firth, file = paste0(out_path, "reg_data/new_setting_labels_firth.RData"))
# save(new_setting_labels_firth, file = paste0(out_path, "reg_data/new_setting_labels_firth_g_2001.RData"))
save(new_setting_labels_firth, file = paste0(out_path, "reg_data/new_setting_labels_firth_no_year_month.RData"))



#### Missed opportunities All (new setting labels - remove sep)-------------------------
library(detectseparation)

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd_2 %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations_2 %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  test <- glm(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult
              # + year + month
              ,
              family = "binomial", data=reg_data, method = "detect_separation")
  
  if(!test$outcome){
    fit <- glm(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult 
               # + year + month
               ,
               family = "binomial", data=reg_data)
    
    out <- tibble(term = names(coef(fit)),
                  estimate = coef(fit))
  } else {
    out <- tibble(term = "quasi_sep",
                  estimate = NA)
  }
  return(out)
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

new_setting_labels_remove <- bind_rows(miss_opp_res) %>% 
  filter(term != "quasi_sep") %>% 
group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

count_quasi_sep <-  bind_rows(miss_opp_res) %>% 
  filter(term == "quasi_sep") %>% 
  count() %>% .$n

new_setting_labels_remove <- bind_rows(new_setting_labels_remove,
                                            tibble(term = "sim_trials_removed",
                                                   est = count_quasi_sep,
                                                   low = NA,
                                                   high = NA))

# save(new_setting_labels_remove, file = paste0(out_path, "reg_data/new_setting_labels_remove.RData"))
# save(new_setting_labels_remove, file = paste0(out_path, "reg_data/new_setting_labels_remove_g_2001.RData"))
save(new_setting_labels_remove, file = paste0(out_path, "reg_data/new_setting_labels_remove_no_year_month.RData"))


################################################################################

#### Missed opportunities All (main effect for setting)-------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  fit <- glm(miss~ inpatient + ed + outpatient + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
             # + year + month
             , 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)

gc()

main_effects_setting_labels <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# table4 <- table4 %>% 
#   mutate(nice = paste0(trimws(format(round(est, 2), nsmall = 2)),
#                        " (", 
#                        trimws(format(round(low, 2), nsmall = 2)),
#                        "-",
#                        trimws(format(round(high, 2), nsmall = 2)),
#                        ")"))
# rm(cluster)
# write_csv(table4, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table4.csv"))
# 

# save(main_effects_setting_labels, file = paste0(out_path, "reg_data/main_effects_setting_labels.RData"))
# save(main_effects_setting_labels, file = paste0(out_path, "reg_data/main_effects_setting_labels_g_2001.RData"))
save(main_effects_setting_labels, file = paste0(out_path, "reg_data/main_effects_setting_labels_no_year_month.RData"))


#### Missed opportunities All (main effects setting labels - Firth regression)-------------------------
library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
   
  # fit <- glm(miss~ inpatient + ed + outpatient + obs_stay + age_cat + female + source + weekend + abx + ID_consult + year + month, 
  #            family = "binomial", data=reg_data)
  
  fit <- logistf(miss~ inpatient + ed + outpatient + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
                 # + year + month
                 ,
                 family = "binomial", data=reg_data,
                 control = logistf.control( maxit = 1000),
                 plcontrol = logistpl.control( maxit = 1000))
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  # broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

main_effects_setting_labels_firth <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# save(main_effects_setting_labels_firth, file = paste0(out_path, "reg_data/main_effects_setting_labels_firth.RData"))
# save(main_effects_setting_labels_firth, file = paste0(out_path, "reg_data/main_effects_setting_labels_firth_g_2001.RData"))
save(main_effects_setting_labels_firth, file = paste0(out_path, "reg_data/main_effects_setting_labels_firth_no_year_month.RData"))



#### Missed opportunities All (main effects setting labels - remove sep)-------------------------
library(detectseparation)

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
    # filter(year > 2001) %>% 
    # mutate(year = factor(year, 
    #                      levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
  
  test <- glm(miss~ inpatient + ed + outpatient + obs_stay + age_cat + female + source + weekend + abx + ID_consult 
              # + year + month
              ,
              family = "binomial", data=reg_data, method = "detect_separation")
  
  if(!test$outcome){
    fit <- glm(miss~ inpatient + ed + outpatient + obs_stay + age_cat + female + source + weekend + abx + ID_consult
               # + year + month
               ,
               family = "binomial", data=reg_data)
    
    out <- tibble(term = names(coef(fit)),
                  estimate = coef(fit))
  } else {
    out <- tibble(term = "quasi_sep",
                  estimate = NA)
  }
  return(out)
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

main_effects_setting_labels_remove <- bind_rows(miss_opp_res) %>% 
  filter(term != "quasi_sep") %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

count_quasi_sep <- bind_rows(miss_opp_res) %>% 
  filter(term == "quasi_sep") %>% 
  count() %>% .$n

main_effects_setting_labels_remove <- bind_rows(main_effects_setting_labels_remove,
                                            tibble(term = "sim_trials_removed",
                                                   est = count_quasi_sep,
                                                   low = NA,
                                                   high = NA))

# save(main_effects_setting_labels_remove, file = paste0(out_path, "reg_data/main_effects_setting_labels_remove.RData"))
# save(main_effects_setting_labels_remove, file = paste0(out_path, "reg_data/main_effects_setting_labels_remove_g_2001.RData"))
save(main_effects_setting_labels_remove, file = paste0(out_path, "reg_data/main_effects_setting_labels_remove_no_year_month.RData"))


#### Missed opportunities All (only main effect for inpatient > 2001) -------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>%
                          mutate(miss=TRUE),
                        index_locations %>%
                          mutate(miss=FALSE)) %>%
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    # mutate(year = factor(year,
    #                      levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
    filter(year > 2001) %>%
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
  
  fit <- glm(miss~inpatient + age_cat + female + source + weekend + abx + opioid + ID_consult + year + month, 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)

gc()

inpatient_ind_only <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


save(inpatient_ind_only, file = paste0(out_path, "reg_data/inpatient_ind_only_g_2001.RData"))



#### Missed opportunities All (new setting labels - Firth regression > 2001)-------------------------
library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd_2 %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations_2 %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    # mutate(year = factor(year,
    #                      levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
    filter(year > 2001) %>%
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
  
  fit <- logistf(miss~setting_label + age_cat + female + source + weekend + abx + opioid + ID_consult  + year + month
                 ,
                 family = "binomial", data=reg_data,
                 control = logistf.control( maxit = 1000),
                 plcontrol = logistpl.control( maxit = 1000))
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  # broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

all_but_setting_vars <- bind_rows(miss_opp_res) %>% 
  filter(!grepl("setting_label",term)) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))


outpatient_data <- parallel::mclapply(miss_opp_res,
                                 function(x){x %>% 
                                     filter(grepl("setting_label",term)) %>%
                                     mutate(term = stringr::str_remove_all(term, "setting_label")) %>% 
                                     bind_rows(tibble(term = "Outpatient Only", estimate = 0))}, 
                                 mc.cores = 60)

outpatient_res <- bind_rows(outpatient_data) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

summary_setting <- function(setting){
  data <- parallel::mclapply(miss_opp_res,
                             function(x){x %>% 
                                 filter(grepl("setting_label",term)) %>%
                                 mutate(term = stringr::str_remove_all(term, "setting_label")) %>% 
                                 bind_rows(tibble(term = "Outpatient Only", estimate = 0)) %>% 
                                 mutate(ref = x %>% 
                                          filter(grepl("setting_label",term)) %>%
                                          mutate(term = stringr::str_remove_all(term, "setting_label")) %>% 
                                          filter(term == setting) %>% .$estimate) %>% 
                                 mutate(estimate = estimate - ref) %>% 
                                 select(-ref)}, 
                             mc.cores = 60)
  
  res <- bind_rows(data) %>% 
    group_by(term) %>% 
    summarise(est = median(exp(estimate), na.rm = T),
              low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
              high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
}

ed_res <- summary_setting("ED Only")
inpatient_res <- summary_setting("Inpatient Only")
obs_res <- summary_setting("Observational Stay Only")

save(miss_opp_res, file = paste0(out_path, "reg_data/new_setting_labels_firth_g_2001_w_diff_ref_setting_data.RData"))
save(all_but_setting_vars,
     outpatient_res,
     ed_res,
     inpatient_res,
     obs_res,
     file = paste0(out_path, "reg_data/new_setting_labels_firth_g_2001_w_diff_ref_setting_summary_res.RData"))


## TESTING Setting ref = "ED"

library(logistf)
get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- sim_res_ssd_2 %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% 
                          mutate(miss=TRUE),
                        index_locations_2 %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    # mutate(year = factor(year,
    #                      levels = paste0("20", str_pad(1:21, width = 2, side = "left", pad = "0"))),
    #        month = factor(month, levels = 1:12))
    filter(year > 2001) %>%
    mutate(year = factor(year,
                         levels = paste0("20", str_pad(2:21, width = 2, side = "left", pad = "0"))),
           month = factor(month, levels = 1:12))
  
  # fit <- glm(miss~setting_label + age_cat + female + source + weekend + abx + ID_consult + year + month, 
  #            family = "binomial", data=reg_data)
  
  reg_data$setting_label <- relevel(reg_data$setting_label, ref = "ED Only")
  fit <- logistf(miss~ setting_label + age_cat + female + source + weekend + abx + ID_consult  + year + month
                 ,
                 family = "binomial", data=reg_data,
                 control = logistf.control( maxit = 1000),
                 plcontrol = logistpl.control( maxit = 1000))
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  # broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(sim_res_ssd$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = 60)

gc()

new_setting_labels_firth_ED_ref <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

save(new_setting_labels_firth_ED_ref, file = paste0(out_path, "reg_data/new_setting_labels_firth_g_2001_w_ED_ref_setting.RData"))

