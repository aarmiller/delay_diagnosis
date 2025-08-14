
# This script provides an example for generating the index visit information using
# the updated system.
#
# There is one dataset generated from this script:
# 1) index_dx_visits.RData - the final index visits that are used 


library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "blasto"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to blasto DB
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

## Stratify by inpatient and outpatient ----------------------------------------

index_cases_org <- index_cases

tmp <- codeBuildr::load_disease_codes("blasto",return_tibble = TRUE)

# collect blasto diagnoses
blasto_dx <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(tmp, copy = TRUE) %>% 
  collect()

# ids with sufficient enrollment
include_ids <- db %>% 
  tbl("index_dx_dates") %>% 
  collect() %>% 
  inner_join(index_cases_org, by = "patient_id")

inpatient_cases <- blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date)) %>% 
  inner_join(blasto_dx) %>% 
  filter(date>=index_dx & date<=(index_dx+7)) %>% 
  filter(inpatient==1) %>% 
  distinct(patient_id) %>% 
  inner_join(select(include_ids,patient_id))

indeterminite_cases <- blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date)) %>% 
  inner_join(blasto_dx) %>% 
  filter(date>index_dx+7 & date<=(index_dx+30)) %>% 
  filter(inpatient==1) %>% 
  distinct(patient_id) %>% 
  anti_join(inpatient_cases) %>% 
  inner_join(select(include_ids,patient_id))

outpatient_cases <- select(include_ids,patient_id) %>% 
  anti_join(inpatient_cases) %>% 
  anti_join(indeterminite_cases)

### Build the index_cases for the blasto_inpatient cohort
cond_name <- "blasto_inpatient"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_org %>% 
  inner_join(inpatient_cases, by = "patient_id")

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))

### Build the index_cases for the blasto_outpatient cohort
cond_name <- "blasto_outpatient"
delay_params <- final_delay_params[[cond_name]]

if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

index_cases <- index_cases_org %>% 
  inner_join(outpatient_cases, by = "patient_id")

save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))

# ## Get location info and build nested cohorts ----------------------------------
# 
# source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
# load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024
# 
# enroll_collapsed_temp <- gather_collapse_enrollment(enrolid_list = index_cases %>% distinct(patient_id) %>% .$patient_id,
#                                                     vars = "egeoloc",
#                                                     db_path =  paste0(delay_params$small_db_path,"blasto.db"),
#                                                     num_cores=10,
#                                                     collect_tab = collect_table(year = 1:22))
# 
# #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322071/
# enroll_collapsed_temp2 <- enroll_collapsed_temp %>% 
#   inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
#   mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
#   inner_join(select(index_cases,patient_id,index_date)) %>% 
#   filter(index_date<=dtend & index_date>=dtstart) %>% 
#   distinct() 
# # only 3,567 out of 3,825 have location information
# 
# # top2_high_inc_state_baddley baddely states with incidence >= 1.38
# 
# location_ind <- enroll_collapsed_temp2 %>%  
#   mutate(top2_high_inc_state_baddley = ifelse(state_abb %in% c("ND", "MN",
#                                                                "WI", "IL",
#                                                                "TN", "AL",
#                                                                "MS", "LA"), 1L, ifelse(is.na(state_abb), NA, 0))) %>% 
#   select(patient_id, location:state_abb, top2_high_inc_state_baddley) %>% 
#   distinct()
# # location_ind %>% filter(top2_high_inc_state_baddley == 1) %>% distinct(state_abb) #check coding
# # 3,387 of the  3,567 with location info have non missing state
# 
# index_cases_w_loc <- index_cases %>%
#   left_join(location_ind %>% select(patient_id, top2_high_inc_state_baddley)) 
# 
# ### Build the index_cases for the blasto_top2_baddley cohort
# cond_name <- "blasto_top2_baddley"
# delay_params <- final_delay_params[[cond_name]]
# 
# if (!dir.exists(delay_params$out_path)) {
#   dir.create(delay_params$out_path)
# }
# 
# index_cases <- index_cases_w_loc %>%
#   filter(top2_high_inc_state_baddley == 1) %>% 
#   select(-top2_high_inc_state_baddley)
#   
# save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))
# 
# ### Build the index_cases for the blasto_not_top2_baddley cohort
# cond_name <- "blasto_not_top2_baddley"
# delay_params <- final_delay_params[[cond_name]]
# 
# if (!dir.exists(delay_params$out_path)) {
#   dir.create(delay_params$out_path)
# }
# 
# index_cases <- index_cases_w_loc %>%
#   filter(top2_high_inc_state_baddley == 0) %>% 
#   select(-top2_high_inc_state_baddley)
# 
# save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))
# 
# ### Build the index_cases for the blasto_NA_state cohort
# cond_name <- "blasto_NA_state"
# delay_params <- final_delay_params[[cond_name]]
# 
# if (!dir.exists(delay_params$out_path)) {
#   dir.create(delay_params$out_path)
# }
# 
# index_cases <- index_cases_w_loc %>%
#   filter(is.na(top2_high_inc_state_baddley)) %>% 
#   select(-top2_high_inc_state_baddley)
# 
# save(index_cases,  file = paste0(delay_params$out, "index_cases.RData"))
