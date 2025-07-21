
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "sepsis_kaiser_cp_14"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params_kaiser.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 14 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Connect to db
con <-  src_sqlite(list.files(delay_params$small_db_path, full.names = T))

### Load index cases
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
            "index_cases.RData"))

num_cores <- 40
cp_selected <- delay_params$cp -1 # minus 1 as the risk factors are for within delay window

#update demo1 and demo2
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
# demo2 <- demo2 %>% select(-index_date) %>% 
#   inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
#                select(patient_id, index_date), by = "patient_id")
# 
# rural_ids <- rural_visits %>% inner_join(index_cases %>% distinct(patient_id)) %>% distinct(patient_id)


# update time map
load(paste0(delay_base_path,"delay_tm.RData"))
tm <- tm %>%   
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>%
  mutate(days_since_index = days_since_index - shift) %>%
  select(-shift) %>%
  filter(days_since_index<=0) %>%
  select(patient_id, days_since_index, outpatient, ed, inpatient, other) %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0)) # subset to only visits from AV, ED, IP, or IS

# update all_dx_visits
load(paste0(delay_base_path,"all_dx_visits.RData"))
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

all_dx_visits <- all_dx_visits %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index")) # visit days from other encounter types only removed

# update sim_obs
sim_obs <- sim_obs %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<0) %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index"))

problem_patient_ids <- tm %>% filter(days_since_index == 0) %>% distinct(patient_id) %>%
  anti_join(tm %>% distinct(patient_id),.)


#Load in the setting types
load(paste0(delay_base_path,"visit_info.RData"))

# # update tm_stdplac
# tm_stdprov <- tm_stdprov %>% 
#   inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
#   mutate(days_since_index = days_since_index - shift) %>% 
#   select(-shift) %>%
#   filter(days_since_index<=0)

load(paste0(sim_out_path,"/sim_res_ssd.RData"))
sim_res_ssd <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  select(trial, boot_trial, res) %>% 
  unnest(res) %>% 
  select(-boot_id)

load(paste0(sim_out_path,"sim_obs_reduced.RData"))
sim_obs_reduced <- sim_obs_reduced %>% 
  mutate(days_since_index = -period) %>% 
  select(-period)

sim_res_sim_obs <- sim_res_ssd %>% 
  inner_join(sim_obs_reduced, by = "obs") 

load(paste0(sim_out_path,"boot_data.RData"))
boot_data <- boot_data %>% select(boot_trial, boot_sample)

n_trials <- nrow(distinct(sim_res_ssd,trial))

# Demo data
reg_demo <- demo1 %>% 
  rename(old_sex = sex) %>% 
  inner_join(tribble(~old_sex,~sex, 
                     "F" ,  "Female", 
                     "M" ,  "Male" ,
                     "X" ,  "Neither Male Nor Female" ,
                     "O" ,  "Other",
                     "U" ,  "Unknown/uncertain/missing"
  )) %>% 
  rename(old_race = race) %>% 
  inner_join(tribble(~old_race,~race, 
                     "HP","Native Hawaiian / Pacific Islander", 
                     "IN", "American Indian / Alaskan Native",
                     "AS", "Asian",
                     "BA", "Black or African American",
                     "WH", "White",
                     "MU", "Multiple races with particular unknown",
                     "OT", "Other, values that do not fit well in any",
                     "other value ", "Other Value",
                     "UN", "Unknown or Not Reported")) %>% 
  mutate(race = fct_relevel(race,"White")) %>% 
  mutate(sex = fct_relevel(sex,"Male")) %>% 
  mutate( age = index_year-dobyr)

# reg_demo <- reg_demo %>% 
#   left_join(demo2 %>% 
#               filter(index_date<=dtend & index_date>=dtstart) %>% 
#               mutate(msa_new = msa %in% c("0","")) %>% 
#               mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
#               mutate(source = as.factor(source)) %>% 
#               distinct(patient_id,source,msa=msa_new),
#             by = "patient_id")

age_cats <- c(-1,17,34,44,54,64,130)

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,sex,age_cat,age,year,month,race) 

tm <- tm %>% rename(admdate = svcdate)

# ### Prepare rx_visits ----------------------------------------------------
# rx_visits <- con %>% tbl("all_rx_visits") %>% 
#   filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
#   collect()
# 
# rx_visits <- rx_visits %>% 
#   inner_join(select(demo1,patient_id,index_date)) %>% 
#   filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))
# 
# gc()
# 
# ### Prepare Drug Indicators ----------------------------------------------------
# red_book <- codeBuildr::load_redbook(redbook_path = "/Shared/Statepi_Marketscan/databases/Truven/redbook.csv")
# 
# tmp1 <- red_book %>% filter(str_detect(tolower(THRGRDS),"immuno")) %>% .$NDCNUM
# tmp2 <- red_book %>% filter(str_detect(tolower(THRDTDS),"prednisone")) %>% .$NDCNUM
# 
# immuno_rx_ids <- rx_visits %>% 
#   filter(ndcnum %in% tmp1) %>% 
#   filter(date<index_date-(cp_selected) & date>=(index_date-(delay_params$upper_bound)) ) %>%
#   # filter(date<index_date-14 & (date+daysupp)>=(index_date-30)) %>% 
#   distinct(patient_id) %>% 
#   mutate(immunosuppressant=1L)
# 
# pred_rx_ids <- rx_visits %>% 
#   filter(ndcnum %in% tmp2) %>% 
#   filter(date<index_date-(cp_selected-1) & date>=(index_date-(delay_params$upper_bound)) ) %>% 
#   # filter(date<index_date-14 & (date+daysupp)>=(index_date-30)) %>% 
#   distinct(patient_id) %>% 
#   mutate(prednisone=1L)
# 
# reg_demo <- reg_demo %>% 
#   left_join(immuno_rx_ids) %>% 
#   left_join(pred_rx_ids) %>% 
#   mutate_at(vars(immunosuppressant,prednisone),~replace_na(.,0L))
# 
# # combine immunosuppressant and prednisone
# reg_demo <- reg_demo %>% 
#   mutate(immuno_comb = as.integer((immunosuppressant==1 | prednisone==1))) 
# 
# ## Prepare antibiotic indicators -----------------------------------------------
# load("/Shared/AML/truven_mind_projects/antibiotic_risk_categories/antibiotics_groupings_new.RData")
# 
# abx_codes <- antibiotic_ndc_groups_new
# abx_rx_vis <- rx_visits %>% 
#   filter(ndcnum %in% abx_codes$ndcnum) 
# 
# reg_demo <- reg_demo %>% 
#   left_join(abx_rx_vis %>% 
#               filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
#               distinct(patient_id) %>% 
#               mutate(abx_window=1L)) %>% 
#   mutate_at(vars(abx_window),~replace_na(.,0L))
# 
# ## Prepare inhaler indicators --------------------------------------------------
# 
# tmp <- read_csv("/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/data/all_inhaler_codes.csv")
# 
# inhaler_rx_vis <- rx_visits %>% 
#   filter(ndcnum %in% tmp$NDCNUM) 
# 
# reg_demo <- reg_demo %>% 
#   left_join(inhaler_rx_vis %>% 
#               filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
#               distinct(patient_id) %>% 
#               mutate(inhaler_window=1L)) %>% 
#   mutate_at(vars(inhaler_window),~replace_na(.,0L))
# 
# ## Add abx and ID consult to time map
# abx_rx_vis <- abx_rx_vis %>% distinct(patient_id, admdate = date)
# inhaler_rx_vis <- inhaler_rx_vis %>% distinct(patient_id, admdate = date)
# 
# tm <- tm %>% rename(admdate = svcdate) %>% wil not need to rename if did above in line 146
#  left_join(abx_rx_vis %>% mutate(abx = 1L), by = c("patient_id", "admdate")) %>% 
#   left_join(inhaler_rx_vis %>% mutate(inhaler = 1L), by = c("patient_id", "admdate")) %>% 
#   mutate_at(vars(abx, inhaler),~replace_na(.,0L)) 
# 
# ## Prepare opioid codes --------------------------------------------------------
# tmp <- codeBuildr::load_rx_codes("opioids")
# 
# tmp <- unique(unlist(tmp,use.names = F))
# 
# opioid_rx_vis <- rx_visits %>% 
#   filter(ndcnum %in% tmp) 
# 
# reg_demo <- reg_demo %>% 
#   left_join(opioid_rx_vis %>% 
#               filter(date<index_date & (date)>=(index_date-cp_selected)) %>% 
#               distinct(patient_id) %>% 
#               mutate(opioid_window=1L)) %>% 
#   mutate_at(vars(opioid_window),~replace_na(.,0L))
# 
# opioid_rx_vis <- opioid_rx_vis %>% distinct(patient_id, admdate = date)
# 
# tm <- tm %>% 
#   left_join(opioid_rx_vis %>% mutate(opioid = 1L), by = c("patient_id", "admdate")) %>% 
#   mutate_at(vars(opioid),~replace_na(.,0L)) 
# 
# ## Prepare any indicator -------------------------------------------------------
# prior_vis_ind <- tm %>% 
#   inner_join(index_cases) %>% 
#   filter(admdate<index_date & admdate>=(index_date-cp_selected)) %>% 
#   distinct(patient_id) %>% 
#   mutate(prior_vis=1L)
# 
# reg_demo <- reg_demo %>% 
#   left_join(prior_vis_ind) %>% 
#   mutate_at(vars(prior_vis),~replace_na(.,0L))
# 
# ## Prepare num before indicator ------------------------------------------------
# num_vis_before_window_ind <- tm %>% 
#   inner_join(index_cases) %>% 
#   filter(admdate<(index_date-cp_selected) & admdate>=(index_date-(delay_params$upper_bound)) ) %>% 
#   distinct(patient_id,admdate) %>% 
#   count(patient_id,name = "num_vis_before")
# 
# reg_demo <- reg_demo %>% 
#   left_join(num_vis_before_window_ind) %>% 
#   mutate_at(vars(num_vis_before),~replace_na(.,0L))

## Prepare fever indicator ------------------------------------------------
cond_temp <- "sepsis"
vital_signs_info <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", cond_temp, "/", cond_temp, "_vital_signs_12sep23final.sas7bdat"))

vital_signs_info <- vital_signs_info %>% 
  rename(admdate = measurement_date,
         patient_id = STUDYID) %>% 
  mutate(admdate = as.integer(admdate) )%>% 
  inner_join(index_cases %>% select(patient_id, index_date)) %>% 
  mutate(days_since_index = admdate - index_date) %>% 
  filter(days_since_index <= 0 & days_since_index >= -delay_params$upper_bound) 

for (x in c(100, 100.4, 101.3)) {
  temp <- vital_signs_info %>% 
    filter(!is.na(temp)) %>%
    filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
    filter(temp >= x) %>% 
    distinct(patient_id, admdate) 
  
  # add or fever dx 
  
  tm <- tm %>% 
    left_join(temp %>% mutate(ind = 1L), by = c("patient_id", "admdate")) %>% 
    mutate_at(vars(ind),~replace_na(.,0L)) 
  
  names(tm)[names(tm) == "ind"] <- paste0("fever_", x)
}


# vital_signs_info %>% 
#   filter(!is.na(temp)) %>%
#   filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
#   filter(days_since_index <0 & days_since_index>=-cp_selected) %>% 
#   distinct(patient_id)
# 
# vital_signs_info %>% 
#   filter(!is.na(temp)) %>%
#   filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
#   filter(days_since_index==0) %>% 
#   distinct(patient_id)
# 
# vital_signs_info %>% 
#   filter(!is.na(temp)) %>%
#   filter(temp >= 100) %>% 
#   filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
#   filter(days_since_index <0 & days_since_index>=-cp_selected) %>% 
#   distinct(patient_id)
# 
# vital_signs_info %>% 
#   filter(!is.na(temp)) %>%
#   filter(temp >= 100) %>% 
#   filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
#   filter(days_since_index==0) %>% 
#   distinct(patient_id)

############################
#### Prepare Visit Info ####
############################

# total number of combinations: \sum_i=1^k choose(n, k) + 1 for no cat
# so here it is sum(choose(4, 1:4)) +1

setting_labels <- expand.grid(outpatient = c("Outpatient", NA),
                    inpatient = c("Inpatient", NA),
                    ed = c("ED", NA),
                    other = c("Other", NA)) %>% 
  unite(., col = setting_label, sep = " + ", na.rm = T, remove =F) %>% 
  mutate(across(outpatient:other, ~ifelse(is.na(.), F, T))) %>% 
  mutate(setting_label = ifelse(setting_label == "", "Not any",
                                ifelse(setting_label == "Inpatient", "Inpatient Only",
                                       ifelse(setting_label == "Outpatient", "Outpatient Only",
                                              ifelse(setting_label == "ED", "ED Only",
                                                     ifelse(setting_label == "Other", "Other Only",
                                                            setting_label))))))
setting_labels <- setting_labels[order(rowSums(setting_labels[,2:5])), ]
setting_labels <- setting_labels %>% mutate(setting_label = factor(setting_label, 
                                 levels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")]),
                                 labels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")])))

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(outpatient==1 | ed==1  | inpatient==1) %>%
  # filter(outpatient==1 | ed==1 | other==1 | inpatient==1) %>%
  # distinct(patient_id,outpatient,ed,other,inpatient,admdate,abx, opioid, inhaler, fever_100, fever_100.4, fever_101.3) %>% 
  distinct(patient_id,outpatient,ed,other,inpatient,admdate, fever_100, fever_100.4, fever_101.3) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  # select(patient_id,outpatient,ed,other,inpatient,weekend,dow, year, month, abx, opioid, inhaler, fever_100, fever_100.4, fever_101.3) %>% 
  select(patient_id,outpatient,ed,other,inpatient,weekend,dow, year, month, fever_100, fever_100.4, fever_101.3) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1  | inpatient==1) %>%
  # filter(outpatient==1 | ed==1 | other==1 | inpatient==1) %>%
  # distinct(obs,patient_id,outpatient,ed,other,inpatient,admdate,abx, opioid, inhaler, fever_100, fever_100.4, fever_101.3) %>%
  distinct(obs,patient_id,outpatient,ed,other,inpatient,admdate, fever_100, fever_100.4, fever_101.3) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  # select(obs,patient_id,outpatient,ed,other,inpatient,weekend,dow, year, month, abx, opioid, inhaler, fever_100, fever_100.4, fever_101.3) %>% 
  select(obs,patient_id,outpatient,ed,other,inpatient,weekend,dow, year, month, fever_100, fever_100.4, fever_101.3) %>%
  inner_join(setting_labels) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd <- sim_res_ssd %>% 
  inner_join(obs_locations)


full_reg_data <- nest(sim_res_ssd, .by = c("trial", "boot_trial"))
full_reg_data <- full_reg_data %>% 
  inner_join(boot_data %>% 
               mutate(boot_sample = map(boot_sample, ~inner_join(., index_locations))))

# all.equal(full_reg_data %>% slice(100) %>% select(boot_sample) %>% unnest(boot_sample),
#           full_reg_data %>% slice(1000) %>% select(boot_sample) %>% unnest(boot_sample))

print("dataset ready")

###########################
#### Regression Models ####
###########################

#### Missed opportunities All --------------------------------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>%
                          unnest(data) %>%
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>%
                          unnest(boot_sample) %>%
                          mutate(miss=FALSE)) %>%
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>%
    select(miss, age_cat, sex, weekend, year, month,
            fever_100) %>%
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  fit <- glm(miss~., 
             family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                            function(x){get_miss_res(x)}, 
                                            mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T))

# rm(cluster)
gc()

# 
# #### Missed opportunities Medicaid ---------------------------------------------
# 
# get_miss_res_med <- function(trial_val){
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo, by = "patient_id") %>% 
#     filter(source == "medicaid")
#   
#   
#   fit <- glm(miss~setting_label + age_cat + female + rural + race + weekend + abx + opioid + ID_consult,
#              family = "binomial", data=reg_data)
#   
#   # fit.penalized <- logistf(miss~setting_label + age_cat + female + rural + race + weekend + abx ,
#   #                          family = binomial, data=reg_data,
#   #                          control = logistf.control( maxit = 1000),
#   #                          plcontrol = logistpl.control( maxit = 1000)) 
#   # 
#   # tibble(term = names(coef(fit.penalized)),
#   #        estimate = coef(fit.penalized))
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# # get_miss_res_med(10)
# 
# # cluster <- parallel::makeCluster(num_cores)
# # 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("get_miss_res_med","sim_res_ssd","index_locations","reg_demo"),
# #                         envir=environment())
# # 
# # 
# # miss_opp_res_med <- parallel::parLapply(cl = cluster,
# #                                         1:max(sim_res_ssd$trial),
# #                                         function(x){get_miss_res_med(x)})
# # 
# # 
# # parallel::stopCluster(cluster)
# 
# 
# miss_opp_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
#                                             function(x){get_miss_res_med(x)}, 
#                                             mc.cores = num_cores)
# 
# miss_opp_res_med <- bind_rows(miss_opp_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# 
# # rm(cluster)
# gc()
# 
# 
# #### Delay duration All --------------------------------------------------------
# 
# # compute duration by simulation (note this will be used in the next step)
# full_reg_data_dur <- full_reg_data %>% 
#   mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")) %>% 
#                       group_by(patient_id) %>% 
#                       summarise(duration = -min(days_since_index)) %>% 
#                       ungroup()))
# 
# # rm(sim_res_ssd)
# gc()
#   
# # 
# # all.equal(full_reg_data %>% select(data) %>% slice(1) %>% 
# #             unnest(data) %>% inner_join(., sim_obs_reduced,by = c("obs", "patient_id")) %>% 
# #             group_by(patient_id) %>% 
# #             summarise(duration = -min(days_since_index)) %>% 
# #             ungroup(),
# #           full_reg_data %>% select(data) %>% slice(1) %>% 
# #             unnest(data) %>% inner_join(., sim_obs,by = c("obs", "patient_id")) %>% 
# #             group_by(patient_id) %>% 
# #             summarise(duration = -min(days_since_index)) %>% 
# #             ungroup())
# 
# get_dur_res <- function(trial_val){
#   
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L))
#   
#   
#   fit <- glm(duration~age_cat + female + rural + source + abx_window + opioid_window, family = "gaussian", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# # get_dur_res(10)
# 
# # cluster <- parallel::makeCluster(num_cores)
# # 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("sim_res_dur","get_dur_res","reg_demo"),
# #                         envir=environment())
# # 
# # 
# # miss_dur_res <- parallel::parLapply(cl = cluster,
# #                                     1:n_trials,
# #                                     function(x){get_dur_res(x)})
# # 
# # parallel::stopCluster(cluster)
# 
# miss_dur_res <- parallel::mclapply(1:max(full_reg_data_dur$trial),
#                                        function(x){get_dur_res(x)}, 
#                                        mc.cores = num_cores)
# 
# 
# miss_dur_res <- bind_rows(miss_dur_res) %>% 
#   group_by(term) %>% 
#   summarise(est = median(estimate, na.rm = T),
#             low = quantile(estimate, probs = c(0.025), na.rm = T),
#             high = quantile(estimate, probs = c(0.975), na.rm = T))
# 
# 
# gc()
# 
# #### Delay duration Medicaid ---------------------------------------------------
# 
# get_dur_res_med <- function(trial_val){
#   
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L)) %>% 
#     filter(source == "medicaid")
#   
#   
#   fit <- glm(duration~age_cat + female + rural + race + abx_window + opioid_window, family = "gaussian", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# # get_dur_res_med(10)
# 
# # cluster <- parallel::makeCluster(num_cores)
# 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("get_dur_res_med"),
# #                         envir=environment())
# # 
# # 
# # miss_dur_res_med <- parallel::parLapply(cl = cluster,
# #                                         1:n_trials,
# #                                         function(x){get_dur_res_med(x)})
# # 
# # parallel::stopCluster(cluster)
# 
# miss_dur_res_med <- parallel::mclapply(1:max(full_reg_data_dur$trial),
#                                    function(x){get_dur_res_med(x)}, 
#                                    mc.cores = num_cores)
# 
# 
# miss_dur_res_med <- bind_rows(miss_dur_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(estimate, na.rm = T),
#             low = quantile(estimate ,probs = c(0.025), na.rm = T),
#             high = quantile(estimate ,probs = c(0.975), na.rm = T))
# 
# gc()
# 
# #### Delay Patient All ---------------------------------------------------------
# 
# 
# get_delay_pat_res <- function(trial_val){
#   
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L)) %>%  
#     mutate(miss = duration>0)
#   
#   
#   fit <- glm(miss~age_cat + female + rural + source + abx_window + opioid_window, family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# # get_delay_pat_res(10)
# 
# # cluster <- parallel::makeCluster(num_cores)
# 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("get_delay_pat_res"),
# #                         envir=environment())
# # 
# # 
# # miss_delay_pat_res <- parallel::parLapply(cl = cluster,
# #                                           1:n_trials,
# #                                           function(x){get_delay_pat_res(x)})
# 
# # parallel::stopCluster(cluster)
# 
# miss_delay_pat_res <- parallel::mclapply(1:max(full_reg_data_dur$trial),
#                                        function(x){get_delay_pat_res(x)}, 
#                                        mc.cores = num_cores)
# 
# miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# gc()
# 
# #### Delay Patient Medicaid ----------------------------------------------------
# 
# 
# get_delay_pat_res_med <- function(trial_val){
# 
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L)) %>%  
#     mutate(miss = duration>0) %>% 
#     filter(source == "medicaid")
#   
#   
#   fit <- glm(miss~age_cat + female + rural + race + abx_window + opioid_window, family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# # get_delay_pat_res_med(10)
# 
# # cluster <- parallel::makeCluster(num_cores)
# 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("get_delay_pat_res_med"),
# #                         envir=environment())
# # 
# # 
# # miss_delay_pat_res_med <- parallel::parLapply(cl = cluster,
# #                                               1:n_trials,
# #                                               function(x){get_delay_pat_res_med(x)})
# # 
# # 
# # parallel::stopCluster(cluster)
# 
# miss_delay_pat_res_med <- parallel::mclapply(1:max(full_reg_data_dur$trial),
#                                          function(x){get_delay_pat_res_med(x)}, 
#                                          mc.cores = num_cores)
# 
# miss_delay_pat_res_med <- bind_rows(miss_delay_pat_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# 
# gc()
# 
# ### Alternative missed visits --------------------------------------------------
# # 
# # obs_locations2 <- tm %>% 
# #   left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
# #   filter(!is.na(obs)) %>% 
# #   filter(setting_type !=4) %>%
# #   distinct(patient_id,obs,setting_type,admdate) %>% 
# #   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
# #   mutate(dow=weekdays(as_date(admdate))) %>% 
# #   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
# #   select(patient_id,obs,setting,dow,weekend) 
# # 
# # index_locations2 <- tm %>% 
# #   filter(days_since_index==0) %>% 
# #   filter(setting_type !=4) %>%
# #   distinct(patient_id,setting_type,admdate) %>% 
# #   mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
# #   mutate(dow=weekdays(as_date(admdate))) %>% 
# #   mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
# #   select(patient_id,setting,weekend,dow) 
# # 
# # sim_res_ssd <- sim_res_ssd %>% 
# #   select(obs,trial) %>% 
# #   inner_join(obs_locations2)
# # 
# # get_miss_res2 <- function(trial_val){
# #   
# #   # trial_val <- 1
# #   
# #   tmp1 <- sim_res_ssd %>% 
# #     filter(trial==trial_val)
# #   
# #   reg_data <- bind_rows(tmp1 %>% 
# #                           mutate(miss=TRUE),
# #                         index_locations2 %>% 
# #                           mutate(miss=FALSE)) %>% 
# #     inner_join(reg_demo, by = "patient_id")
# #   
# #   
# #   fit <- glm(miss~setting + age_cat + female + rural + source + weekend, family = "binomial", data=reg_data)
# #   
# #   broom::tidy(fit)
# #   
# #   
# # }
# # 
# # get_miss_res2(10)
# # 
# # cluster <- parallel::makeCluster(num_cores)
# # # cluster <- parallel::makeCluster(6)
# # 
# # parallel::clusterCall(cluster, function() library(tidyverse))
# # parallel::clusterExport(cluster,c("sim_res_ssd","get_miss_res2","index_locations2","reg_demo"),
# #                         envir=environment())
# # 
# # 
# # miss_opp_res2 <- parallel::parLapply(cl = cluster,
# #                                     1:max(sim_res_ssd$trial),
# #                                     function(x){get_miss_res2(x)})
# # 
# # 
# # parallel::stopCluster(cluster)
# # 
# # miss_opp_res2 <- bind_rows(miss_opp_res2) %>% 
# #   group_by(term) %>% 
# #   summarise(est = mean(exp(estimate)),
# #             low = quantile(exp(estimate),probs = c(0.025)),
# #             high = quantile(exp(estimate),probs = c(0.975)))
# # 
# # 
# # miss_opp_res2
# # miss_opp_res

########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res)
                             # miss_opp_res_med = miss_opp_res_med,
                             # miss_dur_res = miss_dur_res,
                             # miss_dur_res_med = miss_dur_res_med,
                             # miss_delay_pat_res = miss_delay_pat_res,
                             # miss_delay_pat_res_med = miss_delay_pat_res_med)


save(ssd_miss_risk_models,file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/sepsis_kaiser/risk_model_report.Rmd",
                  params = list(cond = "Sepsis Kaiser (cp = 14)"),
                  output_dir = out_path,
                  output_file = paste0(proj_name, "_risk_model_report_", lubridate::today() %>% format('%m-%d-%Y')))

