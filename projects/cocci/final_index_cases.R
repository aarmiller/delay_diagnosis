
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 


library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "cocci"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to cocci DB
db <- src_sqlite(paste0(delay_params$small_db_path, cond_name, ".db"))

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()

# Set enrollment minimum
enroll_min <- delay_params$upper_bound

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

index_cases <- enrolid_crosswalk %>%
  inner_join(index_dx_dates) %>% 
  select(patient_id,enrolid,index_date) %>% 
  mutate(shift=0L)

save(index_cases, file = paste0(delay_params$out, "index_cases.RData"))

## Get AZ location info and build nested cohorts -------------------------------

load("/Shared/Statepi_Diagnosis/projects/cocci/delay_window_1_91/sim_results/AZ_ind_data.RData")
#comes from Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/cocci/sim_res_stratified_by_AZ

index_cases_w_loc <- index_cases %>% 
  left_join(AZ_ind_data %>% select(patient_id, AZ)) 

index_cases_w_loc %>% count(AZ)
AZ_ind_data %>% count(AZ)

### Build the index_cases for the AZ cohort
cond_name <- "cocci_AZ"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_w_loc %>%
  filter(AZ == 1) %>% 
  select(-AZ)

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))

### Build the index_cases for the non-AZ cohort
cond_name <- "cocci_not_AZ"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_w_loc %>%
  filter(AZ == 0) %>% 
  select(-AZ)

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))
