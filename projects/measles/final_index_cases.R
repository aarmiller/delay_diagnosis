
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 
# 2) index_dx_visits_validated.RData - the final index visits that are used for 
#    secondary/sensitivity analysis.


library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "measles"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to measles DB
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
cond_name <- "measles_validated"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

# Aaron Miller 10/10/2025
# 7:38 AM
# Hey Alan, Sorry I missed this.
# Yes, we should be using 86765 instead. 87798 is also correct, this is the code used for PCR testing of an active infection.
# So the codes we should use are "86765" and "87798". I updated the script on the repo and it looks like this brings up the case count.

# CPT 87798 — Infectious agent detection by nucleic acid (DNA or RNA), not otherwise specified; amplified probe technique, each organism.
# CPT 86765 - Measles (Rubeola) Antibody, IgG test

proc_codes <- c("86765","87798")

validated_ids <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% proc_codes) %>% 
  filter(between(days_since_index,-30,30)) %>% 
  collect() %>% 
  group_by(patient_id) %>% 
  summarise(shift = min(days_since_index)) %>% 
  distinct(patient_id,shift)

validated_ids <- validated_ids %>% 
  mutate(shift = ifelse(shift>0,0,shift))

index_cases <- index_cases %>% 
  select(-shift) %>% 
  inner_join(validated_ids)

# index_cases %>% filter(shift != 0) %>%
#   summarise(n = n(),
#            mean_shift = mean(shift),
#            median_shift = median(shift))

# A tibble: 1 x 3
# n mean_shift median_shift
# <int>      <dbl>        <dbl>
#   1   172      -11.2           -9

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))

