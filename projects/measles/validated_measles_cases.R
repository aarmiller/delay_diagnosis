rm(list = ls())

library(tidyverse)

db <- src_sqlite("/Volumes/AML/truven_extracts/small_dbs/measles/measles.db")

index_dx <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

# current cases
index_dx %>% 
  filter(time_before_index>=180)

# find cases with lab testing within 30 days of index diagnosis
validated_cases <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("86765","87798")) %>% 
  filter(between(days_since_index,-30,30)) %>% 
  distinct(patient_id,date) %>% 
  collect()

# find earliest testing date
validated_cases <- validated_cases %>% 
  group_by(patient_id) %>% 
  summarise(testing_date = min(date))
  

final_cases <- index_dx %>% 
  filter(time_before_index>=180) %>% 
  select(patient_id,index_date) %>% 
  inner_join(validated_cases) %>% 
  mutate(new_index = ifelse(index_date<=testing_date,index_date,testing_date))

