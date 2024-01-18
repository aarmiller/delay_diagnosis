
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

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis.db")


## Build the primary index_dx_visits -------------------------------------------
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()
index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

final_index_dates <- enrolid_crosswalk %>%
  inner_join(index_dx_dates) %>% 
  select(patient_id,enrolid,index_date) %>% 
  mutate(shift=0L)

save(final_index_dates, file = "/Shared/Statepi_Diagnosis/projects/pertussis/final_index_dates.RData")

## Build the secondary index_dx_visits_validated -------------------------------

# Add family ID
enrolid_crosswalk <- enrolid_crosswalk %>% 
  mutate(efamid = enrolid %/% 100)

# Load validated ids
validated_ids <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("87798","86615")) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect() %>% 
  group_by(patient_id) %>% 
  summarise(shift = min(days_since_index)) %>% 
  distinct(patient_id,shift)

validated_ids <- validated_ids %>% 
  mutate(shift = ifelse(shift>0,0,shift))

# Find validated family members
add_ids <- enrolid_crosswalk %>% 
  inner_join(validated_ids) %>% 
  distinct(efamid) %>% 
  inner_join(enrolid_crosswalk) %>% 
  anti_join(validated_ids) %>% 
  distinct(patient_id) %>% 
  mutate(shift=0L)

# combine validated and family
tmp <- bind_rows(validated_ids,
                 add_ids)

final_index_dates_validated <- final_index_dates %>% 
  select(-shift) %>% 
  inner_join(tmp)

save(final_index_dates_validated, file = "/Shared/Statepi_Diagnosis/projects/pertussis/final_index_dates_validated.RData")
