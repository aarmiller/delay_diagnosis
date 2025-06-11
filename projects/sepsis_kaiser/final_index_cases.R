
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 


library(tidyverse)
library(bit64)

# Hi Alan,
# 
# Here is the inclusion/analysis criteria we should use for the Kaiser study:
#   
# For the study population: Include patients if at least X% of their visits from AV, ED, IP, or IS in days 0-14 before diagnosis included a temperature reading. I think our baseline model should probably be ~80%. But I also think we should evaluate how results compare when we use a 90% or 100% window.
# 
# Next, we should remove fever from the SSD list used in the simulation.
# 
# Finally, the indicator in the risk model should be created from both diagnosis codes and temperature readings. If a patient had either a diagnosis code for fever or a temperature recorded that was a fever, you should create an indicator of fever.
# 
# Let me know if you have any questions.
# 
# Thanks,
# Aaron


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

# Find patients where X% of their all visits from AV, ED, IP, or IS between 0-14 days prior have a temp reading and
# Temp readings can only come from AV ED IP or IS encounter types

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

# load all_dx visits
all_dx <- tbl(db,"all_dx_visits") %>% 
  collect() %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift)
# anti_join(all_dx, index_cases) check

# subset to visits that did not come from other encounter type
load(paste0(delay_params$base_path,"delay_results/delay_tm.RData"))
tm <- tm %>%   
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>%
  mutate(days_since_index = days_since_index - shift) %>%
  select(-shift) %>%
  filter(days_since_index<=0) %>%
  select(patient_id, days_since_index, outpatient, ed, inpatient, other) %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0)) # subset to only visits from AV, ED, IP, or IS

all_dx <- all_dx %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index")) # visit days from other encounter types only removed
# all_dx %>% distinct(patient_id, encounter_date, encounter_type) %>%  # check if we remove visits without AV, ED, IP, or IS
#   mutate(ind = 1) %>% 
#   pivot_wider(id_cols = c("patient_id", "encounter_date"),
#               names_from = "encounter_type",
#               values_from = ind) %>% 
#   mutate(across(AV:LO, ~replace_na(., 0L))) %>% 
#   filter(AV+ED+IP+IS ==0)

temp1 <- all_dx %>% 
  filter(days_since_index <= 0 & days_since_index >= -delay_params$cp+1) %>% 
  distinct(patient_id, days_since_index)  
# temp1 %>% distinct(patient_id) number of patients with at least 1 visit form AV, ED, IP, or IS during window [14,0] 
# anti_join(temp1, index_cases) check

temp2 <- temp1  %>%
  left_join(vital_signs_info,
            by = c("patient_id", "days_since_index") ) %>%
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L))

percent_complete <- temp2 %>%
  group_by(patient_id) %>%
  summarise(percent = mean(temp_reading)*100) %>%
  filter(percent>0)

# alternative
# temp2 <- temp1  %>% 
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>% 
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
#   filter(temp_reading == 1) %>% 
#   distinct(patient_id,days_since_index, temp_reading) 
# 
# percent_complete1 <- temp1 %>% 
#   left_join(temp2,
#             by = c("patient_id", "days_since_index")) %>% 
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>% 
#   count(patient_id, temp_reading) %>% 
#   group_by(patient_id) %>% 
#   mutate(percent = n/sum(n)*100) %>% 
#   ungroup() %>% 
#   filter(temp_reading == 1) 

n_patients <- tibble('Percent of SSD visits in window with temp reading AND temp reading on index' = paste0(">=", seq(50, 100, by = 10), "%"),
        n= sapply(seq(50, 100, by = 10), function(x){sum(percent_complete$percent >= x)}))

names(n_patients)[2] <- paste0("Number of patients [", -1*(-delay_params$cp+1),", 0) delay window")

write_csv(n_patients,  file = paste0(delay_params$out, "study_population_counts.csv"))

# Filter to patient that had temp reading on index and 
# 100% of their SSD visits during window also had temp reading

study_pop <- percent_complete %>% 
  filter(percent >= 100) %>% 
  distinct(patient_id)
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

# Find patients where X% of their all visits from AV, ED, IP, or IS between 0-14 days prior have a temp reading and
# Temp readings can only come from AV ED IP or IS encounter types

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

# load all_dx visits
all_dx <- tbl(db,"all_dx_visits") %>% 
  collect() %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift)
# anti_join(all_dx, index_cases) check

# subset to visits that did not come from other encounter type
load(paste0(delay_params$base_path,"delay_results/delay_tm.RData"))
tm <- tm %>%   
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>%
  mutate(days_since_index = days_since_index - shift) %>%
  select(-shift) %>%
  filter(days_since_index<=0) %>%
  select(patient_id, days_since_index, outpatient, ed, inpatient, other) %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0)) # subset to only visits from AV, ED, IP, or IS

all_dx <- all_dx %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index")) # visit days from other encounter types only removed
# all_dx %>% distinct(patient_id, encounter_date, encounter_type) %>%  # check if we remove visits without AV, ED, IP, or IS
#   mutate(ind = 1) %>% 
#   pivot_wider(id_cols = c("patient_id", "encounter_date"),
#               names_from = "encounter_type",
#               values_from = ind) %>% 
#   mutate(across(AV:LO, ~replace_na(., 0L))) %>% 
#   filter(AV+ED+IP+IS ==0)

temp1 <- all_dx %>% 
  filter(days_since_index <= 0 & days_since_index >= -delay_params$cp+1) %>% 
  distinct(patient_id, days_since_index)  
# temp1 %>% distinct(patient_id) number of patients with at least 1 visit form AV, ED, IP, or IS during window [14,0] 
# anti_join(temp1, index_cases) check

temp2 <- temp1  %>%
  left_join(vital_signs_info,
            by = c("patient_id", "days_since_index") ) %>%
  mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L))

percent_complete <- temp2 %>%
  group_by(patient_id) %>%
  summarise(percent = mean(temp_reading)*100) %>%
  filter(percent>0)

# alternative
# temp2 <- temp1  %>% 
#   left_join(vital_signs_info,
#             by = c("patient_id", "days_since_index") ) %>% 
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>%
#   filter(temp_reading == 1) %>% 
#   distinct(patient_id,days_since_index, temp_reading) 
# 
# percent_complete1 <- temp1 %>% 
#   left_join(temp2,
#             by = c("patient_id", "days_since_index")) %>% 
#   mutate(temp_reading = ifelse(is.na(temp_reading), 0L, 1L)) %>% 
#   count(patient_id, temp_reading) %>% 
#   group_by(patient_id) %>% 
#   mutate(percent = n/sum(n)*100) %>% 
#   ungroup() %>% 
#   filter(temp_reading == 1) 

n_patients <- tibble('Percent of SSD visits in window with temp reading AND temp reading on index' = paste0(">=", seq(50, 100, by = 10), "%"),
                     n= sapply(seq(50, 100, by = 10), function(x){sum(percent_complete$percent >= x)}))

names(n_patients)[2] <- paste0("Number of patients [", -1*(-delay_params$cp+1),", 0) delay window")

write_csv(n_patients,  file = paste0(delay_params$out, "study_population_counts.csv"))

# Filter to patient that had temp reading on index and 
# 100% of their SSD visits during window also had temp reading

study_pop <- percent_complete %>% 
  filter(percent >= 100) %>% 
  distinct(patient_id)
# anti_join(study_pop, index_cases) check

index_cases <- index_cases %>% 
  inner_join(study_pop, by = "patient_id")

save(index_cases, file = paste0(delay_params$out, "index_cases.RData"))