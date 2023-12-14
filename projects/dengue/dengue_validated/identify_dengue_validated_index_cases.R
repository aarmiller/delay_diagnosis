
rm(list = ls())
library(tidyverse)

# load final delay param
cond_name <- "dengue_validated"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]

# Connect to dengue DB
db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/dengue/dengue.db")

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()

# Set enrollment minimum
enroll_min <- delay_params$upper_bound

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

#### Identify dengue cases -----------------------------------------------
# Inclusion Procedures
proc_codes <- c("86790", "87449", "87798")

procs <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% proc_codes) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect()

index_cases <- index_dx_dates %>% inner_join(distinct(procs, patient_id), by = "patient_id")
save(index_cases, file = "/Shared/Statepi_Diagnosis/projects/dengue/dengue_validated/index_cases.RData")


