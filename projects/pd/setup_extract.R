

rm(list = ls())
library(tidyverse)
library(bit64)

# load PD data from Jacob and save to AML 
pd_cases <- read_csv("~/Data/MarketScan/truven_extracts/small_dbs/pd/from_jacob/pd_case_def.csv")
pd_cases <- mutate(pd_cases,enrolid = as.integer64(enrolid))
save(pd_cases, file = "/Volumes/AML/truven_extracts/small_dbs/pd/from_jacob/pd_cases.RData")

pd_study_cohort <- read_csv("~/Data/MarketScan/truven_extracts/small_dbs/pd/from_jacob/pd_control_enrolid.csv")
pd_study_cohort <- mutate(pd_study_cohort,enrolid = as.integer64(enrolid))
save(pd_study_cohort, file = "/Volumes/AML/truven_extracts/small_dbs/pd/from_jacob/pd_study_cohort.RData")

# save crosswalk to small_db directory

pd_index_dates <- pd_study_cohort %>% 
  filter(pd==TRUE) %>% 
  select(enrolid,index_date = first_pd_date)

enrolids <- pd_index_dates %>% 
  arrange(enrolid) %>% 
  mutate(medicaid=0L) %>% 
  mutate(patient_id = row_number()) %>% 
  select(enrolid,medicaid,patient_id)

save(enrolids, file = "/Volumes/AML/truven_extracts/small_dbs/pd/pd_enrolid_crosswalk.RData")

## Load Enrollment and save index_dx_date file -----------

## THIS IS DONE ON HPC

# load enrollment
enroll_db <- src_sqlite("/Shared/Statepi_Marketscan/databases/Truven/enrollment_dbs/all_enroll.db")

enroll_periods <- enroll_db %>%
  tbl("ccae_mdcr_collapse_enroll") %>%
  collect()

load("/Shared/AML/truven_extracts/small_dbs/pd/from_jacob/pd_study_cohort.RData")
load("/Shared/AML/truven_extracts/small_dbs/pd/pd_enrolid_crosswalk.RData")

index_dx_date <- pd_study_cohort %>% 
  filter(pd==TRUE) %>% 
  select(enrolid,index_date = first_pd_date) %>% 
  inner_join(enrolids, by = "enrolid") %>% 
  select(enrolid,medicaid,patient_id,index_date)

# comput the continuous amount of time enrolled before
continuous_time_before <- index_dx_date %>%
  left_join(enroll_periods, by = "enrolid") %>%
  filter(dtstart<=index_date & dtend>=index_date) %>%
  mutate(time_before_index=index_date-dtstart) %>%
  select(enrolid, time_before_index)

# compute time before index any
max_time_before <- enroll_periods %>%
  inner_join(distinct(index_dx_date,enrolid)) %>%
  group_by(enrolid) %>%
  filter(dtstart==min(dtstart)) %>%
  ungroup() %>%
  left_join(index_dx_date, by = "enrolid") %>%
  mutate(max_time_before_index=index_date-dtstart) %>%
  select(enrolid, max_time_before_index)

# finalize index_dx_dates
index_dx_date <- index_dx_date %>%
  left_join(continuous_time_before) %>%
  left_join(max_time_before)

save(index_dx_date, file = "/Shared/AML/truven_extracts/small_dbs/pd/pd_index_dx_dates.RData")


qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh pd

