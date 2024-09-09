
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 


library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "sarcoid"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to sarcoid DB
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

save(index_cases, file = paste0(delay_params$out_path, "index_cases.RData"))

index_cases_main <- index_cases

## Build the sarcoid_lung index_dx_visits --------------------------------------

load("/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/sarcoid_lung_patids.RData")

cond_name <- "sarcoid_lung"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_main %>% 
  inner_join(sarcoid_lung, by = "patient_id")

save(index_cases, file = paste0(delay_params$out_path, "index_cases.RData"))

## Build the sarcoid_skin index_dx_visits --------------------------------------

load("/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_skin/sarcoid_skin_patids.RData")

cond_name <- "sarcoid_skin"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_main %>% 
  inner_join(sarcoid_skin, by = "patient_id")

save(index_cases, file = paste0(delay_params$out_path, "index_cases.RData"))
