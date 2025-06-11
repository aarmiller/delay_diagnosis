
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 


library(tidyverse)
library(bit64)

# Get index cases for delay window [14,0) ######################################
# load final delay param
cond_name <- "sepsis_kaiser_cp_14"
load("/Shared/AML/params/final_delay_params_kaiser.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path, recursive = T)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to blasto DB
db <- src_sqlite(list.files(delay_params$small_db_path, full.names = T))

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()

# Set enrollment minimum
enroll_min <- delay_params$upper_bound

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

index_cases <- index_dx_dates %>% 
  select(patient_id,index_date) %>% 
  mutate(shift=0L)

# Find patients where X% of their SSD visits during window have a temp reading and
# and temp reading on index date. Temp readings can only come from 
# AV ED IP or IS encounter types

# Load vital signs information and subset to enrollees with sufficient enrollment time
vital_signs_info <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", str_remove(delay_params$cond_name, "_kaiser"), "/",
                                           str_remove(delay_params$cond_name, "_kaiser"), "_vital_signs_12sep23final.sas7bdat"))

vital_signs_info <- vital_signs_info %>% 
  rename(patient_id = STUDYID) %>% 
  inner_join(index_cases %>% select(patient_id, index_date), by = "patient_id") %>% 
  mutate(days_since_index = as.numeric(measurement_date - index_date)) %>% 
  filter(days_since_index <= 0 & days_since_index >= -delay_params$upper_bound )
  
# Filter to visit days with non-missing temp and encounter_type %in% c("AV", "ED", "IP", "IS")
vital_signs_info <- vital_signs_info %>% 
  filter(!is.na(temp)) %>% 
  filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
  distinct(patient_id, days_since_index) %>% 
  mutate(temp_reading = 1)
# anti_join(vital_signs_info, index_cases) check

# find temp readings on index
temp_index <- vital_signs_info %>% 
  filter(days_since_index == 0) %>% 
  distinct(patient_id)
# anti_join(temp_index, index_cases) check

# load all_dx visits
all_dx <- tbl(db,"all_dx_visits") %>% 
  collect() %>% 
  inner_join(index_cases %>% select(patient_id), by = "patient_id") 
# anti_join(all_dx, index_cases) check

# load ssd codes
ssd_codes <- codeBuildr::load_ssd_codes(delay_params$ssd_name) %>% 
  filter(type %in% c("icd9","icd10")) %>% 
  mutate(dx_ver = ifelse(type=="icd9","09","10")) %>% 
  select(dx = code,dx_ver)

temp1 <- all_dx %>% 
  inner_join(ssd_codes, by = c("dx", "dx_ver")) %>% 
  filter(days_since_index < 0 & days_since_index >= -delay_params$cp+1) %>% 
  distinct(patient_id, days_since_index)  
# temp1 %>% distinct(patient_id) number of patients with at least 1 SSD visit during delay window
# anti_join(temp1, index_cases) check

temp2 <- temp1  %>% 
  left_join(vital_signs_info,
            by = c("patient_id", "days_since_index") ) %>% 
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
  filter(temp_reading == 1) %>% 
  distinct(patient_id,days_since_index, temp_reading) 

percent_complete <- temp1 %>% 
  left_join(temp2,
            by = c("patient_id", "days_since_index")) %>% 
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>% 
  count(patient_id, temp_reading) %>% 
  group_by(patient_id) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  ungroup() %>% 
  filter(temp_reading == 1) %>%
  inner_join(temp_index)

# alternative
# test1 <- temp1  %>%
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>%
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L))
# 
# test1 %>%
#   group_by(patient_id) %>%
#   summarise(percent_w_temp = mean(temp_reading)) %>%
#   filter(percent_w_temp>0)

n_patients <- tibble('Percent of SSD visits in window with temp reading AND temp reading on index' = paste0(">=", seq(50, 100, by = 10), "%"),
        n= sapply(seq(50, 100, by = 10), function(x){sum(percent_complete$percent >= x)}))

names(n_patients)[2] <- paste0("Number of patients [", -1*(-delay_params$cp+1),", 0) delay window")

write_csv(n_patients,  file = paste0(delay_params$out, "study_population_counts.csv"))

# Filter to patient that had temp reading on index and 
# 100% of their SSD visits during window also had temp reading

study_pop <- percent_complete %>% 
  filter(percent >= 100) %>% 
  distinct(patient_id)
       
study_pop_no_ssd <- index_cases %>% 
  anti_join(., temp1 %>% distinct(patient_id), by = "patient_id") %>% # find those patients without an ssd visit during delay window
  inner_join(temp_index, by = "patient_id") %>% # make sure they had temp reading on index
  distinct(patient_id)

# Check
# study_pop_no_ssd %>% inner_join(all_dx %>% 
#                                   inner_join(ssd_codes, by = c("dx", "dx_ver")) %>% 
#                                   filter(days_since_index < 0 & days_since_index >= -delay_params$cp+1))

# study_pop <- bind_rows(study_pop, study_pop_no_ssd) 
# anti_join(study_pop, index_cases) check

index_cases <- index_cases %>% 
  inner_join(study_pop, by = "patient_id")

save(index_cases, file = paste0(delay_params$out, "index_cases.RData"))

# Get index cases for delay window [7,0) ######################################
rm(list=ls())

# load final delay param
cond_name <- "sepsis_kaiser_cp_7"
load("/Shared/AML/params/final_delay_params_kaiser.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path, recursive = T)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to blasto DB
db <- src_sqlite(list.files(delay_params$small_db_path, full.names = T))

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()

# Set enrollment minimum
enroll_min <- delay_params$upper_bound

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

index_cases <- index_dx_dates %>% 
  select(patient_id,index_date) %>% 
  mutate(shift=0L)

# Find patients where X% of their SSD visits during window have a temp reading and
# and temp reading on index date. Temp readings can only come from 
# AV ED IP or IS encounter types

# Load vital signs information and subset to enrollees with sufficient enrollment time
vital_signs_info <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", str_remove(delay_params$cond_name, "_kaiser"), "/",
                                           str_remove(delay_params$cond_name, "_kaiser"), "_vital_signs_12sep23final.sas7bdat"))

vital_signs_info <- vital_signs_info %>% 
  rename(patient_id = STUDYID) %>% 
  inner_join(index_cases %>% select(patient_id, index_date), by = "patient_id") %>% 
  mutate(days_since_index = as.numeric(measurement_date - index_date)) %>% 
  filter(days_since_index <= 0 & days_since_index >= -delay_params$upper_bound )

# Filter to visit days with non-missing temp and encounter_type %in% c("AV", "ED", "IP", "IS")
vital_signs_info <- vital_signs_info %>% 
  filter(!is.na(temp)) %>% 
  filter(encounter_type %in%  c("AV", "ED", "IP", "IS")) %>% 
  distinct(patient_id, days_since_index) %>% 
  mutate(temp_reading = 1)
# anti_join(vital_signs_info, index_cases) check

# find temp readings on index
temp_index <- vital_signs_info %>% 
  filter(days_since_index == 0) %>% 
  distinct(patient_id)
# anti_join(temp_index, index_cases) check

# load all_dx visits
all_dx <- tbl(db,"all_dx_visits") %>% 
  collect() %>% 
  inner_join(index_cases %>% select(patient_id), by = "patient_id") 
# anti_join(all_dx, index_cases) check

# load ssd codes
ssd_codes <- codeBuildr::load_ssd_codes(delay_params$ssd_name) %>% 
  filter(type %in% c("icd9","icd10")) %>% 
  mutate(dx_ver = ifelse(type=="icd9","09","10")) %>% 
  select(dx = code,dx_ver)

temp1 <- all_dx %>% 
  inner_join(ssd_codes, by = c("dx", "dx_ver")) %>% 
  filter(days_since_index < 0 & days_since_index >= -delay_params$cp+1) %>% 
  distinct(patient_id, days_since_index)  
# temp1 %>% distinct(patient_id) number of patients with at least 1 SSD visit during delay window
# anti_join(temp1, index_cases) check

temp2 <- temp1  %>% 
  left_join(vital_signs_info,
            by = c("patient_id", "days_since_index") ) %>% 
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
  filter(temp_reading == 1) %>% 
  distinct(patient_id,days_since_index, temp_reading) 

percent_complete <- temp1 %>% 
  left_join(temp2,
            by = c("patient_id", "days_since_index")) %>% 
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>% 
  count(patient_id, temp_reading) %>% 
  group_by(patient_id) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  ungroup() %>% 
  filter(temp_reading == 1) %>%
  inner_join(temp_index)

# alternative
# test1 <- temp1  %>%
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>%
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L))
# 
# test1 %>%
#   group_by(patient_id) %>%
#   summarise(percent_w_temp = mean(temp_reading)) %>%
#   filter(percent_w_temp>0)

# Compare counts with 14 day window
# temp1_test <- all_dx %>% 
#   inner_join(ssd_codes, by = c("dx", "dx_ver")) %>% 
#   filter(days_since_index < 0 & days_since_index >= -14) %>% 
#   distinct(patient_id, days_since_index)  
# temp1_test %>% distinct(patient_id)
# temp1%>% distinct(patient_id)
# 
# p1 <- temp1  %>%
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>%
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
#   group_by(patient_id) %>%
#   summarise(n = n(),
#             percent_w_temp = mean(temp_reading)) %>%
#   filter(percent_w_temp>=1) %>% 
#   inner_join(temp_index)
# 
# p2 <- temp1_test  %>%
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>%
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
#   group_by(patient_id) %>%
#   summarise(n = n(), 
#             percent_w_temp = mean(temp_reading)) %>%
#   filter(percent_w_temp>=1)%>% 
#   inner_join(temp_index)
# 
# p2 %>% anti_join(p1 %>% distinct(patient_id)) %>% distinct()
# p1 %>% anti_join(p2 %>% distinct(patient_id)) %>% distinct()
# p1 %>% inner_join(p2 %>% distinct(patient_id)) %>% distinct()
# 
# temp1_test %>% 
#   inner_join(p2 %>% anti_join(p1 %>% distinct(patient_id)) %>% distinct(patient_id)) %>% .$days_since_index %>% 
#   range()
# 
# temp1_test %>% 
#   inner_join(p1 %>% anti_join(p2 %>% distinct(patient_id)) %>% distinct(patient_id)) %>% .$days_since_index %>% 
#   range()


n_patients <- tibble('Percent of SSD visits in window with temp reading AND temp reading on index' = paste0(">=", seq(50, 100, by = 10), "%"),
                     n= sapply(seq(50, 100, by = 10), function(x){sum(percent_complete$percent >= x)}))

names(n_patients)[2] <- paste0("Number of patients [", -1*(-delay_params$cp+1),", 0) delay window")

write_csv(n_patients,  file = paste0(delay_params$out, "study_population_counts.csv"))

# Filter to patient that had temp reading on index and 
# 100% of their SSD visits during window also had temp reading

study_pop <- percent_complete %>% 
  filter(percent >= 100) %>% 
  distinct(patient_id)

study_pop_no_ssd <- index_cases %>% 
  anti_join(., temp1 %>% distinct(patient_id), by = "patient_id") %>% # find those patients without an ssd visit during delay window
  inner_join(temp_index, by = "patient_id") %>% # make sure they had temp reading on index
  distinct(patient_id)

# Check
# study_pop_no_ssd %>% inner_join(all_dx %>% 
#                                   inner_join(ssd_codes, by = c("dx", "dx_ver")) %>% 
#                                   filter(days_since_index < 0 & days_since_index >= -delay_params$cp+1))

# study_pop <- bind_rows(study_pop, study_pop_no_ssd)
# anti_join(study_pop, index_cases) check

index_cases <- index_cases %>% 
  inner_join(study_pop, by = "patient_id")

save(index_cases, file = paste0(delay_params$out, "index_cases.RData"))