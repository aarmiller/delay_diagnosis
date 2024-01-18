
rm(list = ls())
library(tidyverse)
library(bit64)


load(paste0("/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/delay_results/all_dx_visits.RData"))

cond_name <- "meningitis_bacterial"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)


# Load index dates
load(paste0(delay_params$base_path,"delay_results/all_dx_visits.RData"))

# pull in the original diagnoses
load(paste0(delay_params$small_db_path,"meningitis_enrolid_crosswalk.RData"))
load("/Shared/AML/truven_extracts/dx/meningitis/meningitis_dx_visits.RData")

all_facility_visits

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
  distinct(enrolid,date) %>% 
  inner_join(enrolids)

grace_period <- 14

tmp_index <- index_dx_dates %>% 
  inner_join(tmp_index_dates) %>% 
  filter(index_date<=date & index_date+grace_period>=date) %>% 
  distinct(patient_id,index_date,time_before_index,max_time_before_index) %>% 
  filter(time_before_index>=delay_params$upper_bound)

# connect to db
db <- src_sqlite(paste0(delay_params$small_db_path,"meningitis.db"))

proc_visits <- db %>% 
  tbl("all_proc_visits") %>% 
  inner_join(select(tmp_index,patient_id),copy = TRUE) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect()

proc_visits %>% count(days_since_index) %>% .$days_since_index


proc_valid_list <- c("62270","70450","70553","70551","70470","70496","8891",
                     "8703","70460","70552","70546","70545")

proc_visits %>% filter(proc %in% proc_valid_list) %>% distinct(patient_id)

proc_visits %>% count(proc) %>% arrange(desc(n)) %>% write_csv(file = "/Shared/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/proc_counts_near_index.csv")
          
tmp_labels <- read_csv("/Shared/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/top_1000_find_a_code.csv")

select(tmp_labels, proc = description, code_type = status)

proc_visits %>% count(proc) %>% arrange(desc(n)) %>% slice(1:1000) %>% 
  left_join(select(tmp_labels, proc = description, code_type = status)) %>% 
  write_csv("/Shared/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/top_1000_labels.csv")

# run local
tmp <- read_csv("/Volumes/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/proc_counts_near_index.csv")

read_csv("~/Downloads/code-status.csv") %>% write_csv("/Volumes/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/top_1000_find_a_code.csv")


### Specific 


#### OLD (from sepsis remove)

# filter to specific index cases
# 1) 180 days of continuous enrollment
# 2) all prior data in the ICD-10 era (180 days after Jan 1 2016)

# index_cases <- index_dx_dates %>% 
#   filter(time_before_index>=delay_params$upper_bound) %>% 
#   filter(index_date>=as.integer(ymd("2016-01-01")))
# 
# save(index_cases, file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/index_cases.RData")

