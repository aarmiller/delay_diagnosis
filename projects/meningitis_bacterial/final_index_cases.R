# This script provides an example for generating the index visit information using
# the updated system for meningitis
#
# There is one dataset generated from this script: index_dx_visits.RData 
# (Note: here we are validating the index cases using procedure codes from the 
# start, therefore only one file is produced):


library(tidyverse)
library(bit64)


grace_period <- 14 # this is the period of time we are allowing to take place between
                   # the index meningitis date (i.e., first time any meningitis code
                   # was recorded and when a bacterial code was recorded)

enroll_before <- 180

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/meningitis/meningitis.db")
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()
index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

#### Find specific bacterial meningitis codes #####

# pull in the original diagnoses
load("/Shared/AML/truven_extracts/dx/meningitis/meningitis_dx_visits.RData")

index_codes <- bind_rows(tibble(dx_ver=9,
                                dx = c("0360","3200","3201","3202","3207",
                                       "32082","32089","3209")),
                         tibble(dx_ver=10,
                                dx = c("A390","G00","G000","G001","G002","G008",
                                       "G009","G01"))) %>% 
  inner_join(codeBuildr::all_icd_labels)

tmp_index_dates <- bind_rows(all_facility_visits %>% 
                               inner_join(index_codes) %>% 
                               distinct(enrolid,date=svcdate),
                             all_inpatient_visits %>% 
                               inner_join(index_codes) %>% 
                               distinct(enrolid,date=admdate),
                             all_outpatient_visits %>% 
                               inner_join(index_codes) %>% 
                               distinct(enrolid,date=svcdate),
                             inpatient_services_visits %>% 
                               inner_join(index_codes) %>% 
                               distinct(enrolid,date=svcdate)) %>% 
  distinct(enrolid,date)

# find cases with bacterial meningitis diagnosis within grace period
tmp_index <- index_dx_dates %>% 
  inner_join(enrolid_crosswalk) %>% 
  inner_join(tmp_index_dates) %>% 
  filter(index_date<=date & index_date+grace_period>=date) %>% 
  distinct(patient_id,index_date,time_before_index,max_time_before_index) %>% 
  filter(time_before_index>=180) %>% 
  select(patient_id,index_date)

#### Validate with procedure codes ####

# connect to db

proc_visits <- db %>% 
  tbl("all_proc_visits") %>% 
  inner_join(select(tmp_index,patient_id),copy = TRUE) %>% 
  filter(between(days_since_index,-grace_period,grace_period)) %>% 
  collect()

proc_valid_list <- c("62270")

proc_validated <- proc_visits %>% 
  filter(proc %in% proc_valid_list) %>% 
  group_by(patient_id) %>% 
  summarise(shift = min(days_since_index)) %>% 
  mutate(shift = ifelse(shift>0,0,shift))

proc_validated

final_index_dates <- index_dx_dates %>% 
  select(patient_id,index_date) %>% 
  inner_join(enrolid_crosswalk) %>% 
  inner_join(proc_validated) %>% 
  select(patient_id,enrolid,index_date,shift)

save(final_index_dates, file = "/Shared/Statepi_Diagnosis/projects/meningitis_bacterial/final_index_dates.RData")



