
# This script provides an example for generating the index visit information using
# the updated system.
#
# There are two datasets generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 
#
# 2) index_dx_visits_validated.RData - the final index visits that are used for 
#    secondary/sensitivity analysis.
#
# Note: In the event that the "validated" cases are to be used for the primary 
# analysis and there are no secondary/sensitivity analysis conducted, the primary
# dataset index_dx_visits.RData should be the only file generated.

library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "dengue"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to dengue DB
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

## Build the secondary index_dx_visits_validated -------------------------------
cond_name <- "dengue_validated"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

proc_codes <- c("86790", "87449", "87798")

validated_ids <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% proc_codes) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect() %>% 
  group_by(patient_id) %>% 
  summarise(shift = min(days_since_index)) %>% 
  distinct(patient_id,shift)

validated_ids <- validated_ids %>% 
  mutate(shift = ifelse(shift>0,0,shift))

index_cases <- index_cases %>% 
  select(-shift) %>% 
  inner_join(validated_ids)

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))

