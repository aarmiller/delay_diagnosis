
library(tidyverse)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis.db")

# Load Index Visits and Crosswalk
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()

index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

# Add family ID
enrolid_crosswalk <- enrolid_crosswalk %>% 
  mutate(efamid = enrolid %/% 100)

# Load validated ids
validated_ids <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("87798","86615")) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect() %>% 
  distinct(patient_id)

# Find validated family members
all_validated_ids <- enrolid_crosswalk %>% 
  inner_join(validated_ids) %>% 
  distinct(efamid) %>% 
  inner_join(enrolid_crosswalk)
  
