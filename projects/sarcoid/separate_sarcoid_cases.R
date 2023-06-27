
library(tidyverse)


# Connect to sarcoid DB
db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/sarcoid/sarcoid.db")

# Collect index dates
index_dx_dates <- tbl(db,"index_dx_dates") %>% collect()

# Set enrollment minimum
enroll_min <- 365*2

# Filter to enrollees with sufficient enrollment time
index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=enroll_min)

#### Identify Sarcoid Lung Cases -----------------------------------------------


# Inclusion Procedures
lung_proc_codes <- c("31628", "31624", "31629", "3324", "31622", "31620", "31623",
                     "31633", "31632", "3323", "3322", "31645", "31625", "88172",
                     "3324", "10022", "76942", "32405", "38510", "38505", "3328",
                     "38525", "3326", "4011")

lung_procs <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% lung_proc_codes) %>% 
  filter(between(days_since_index,-30,0)) %>% 
  collect()

lung_valid_cases <- index_dx_dates %>% inner_join(distinct(lung_procs,patient_id))


#### Identify Sarcoid Skin Cases -----------------------------------------------

# Inclusion Procedures
skin_proc_codes <- c("11100","11101")

skin_procs <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% skin_proc_codes) %>% 
  filter(between(days_since_index,-30,0)) %>% 
  collect()

skin_valid_cases <- index_dx_dates %>% inner_join(distinct(skin_procs,patient_id))

#### Save Results --------------------------------------------------------------

index_cases <- lung_valid_cases

save(index_cases, file = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/index_cases.RData")

index_cases <- skin_valid_cases

save(index_cases, file = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_skin/index_cases.RData")

#### Analyze both cases --------------------------------------------------------

n_enroll <- index_dx_dates %>% 
  distinct(patient_id) %>% 
  nrow()

n_lung <- lung_valid_cases %>% 
  distinct(patient_id) %>% 
  nrow()

n_skin <- skin_valid_cases %>% 
  distinct(patient_id) %>% 
  nrow()

n_overlap <- skin_valid_cases %>% 
  distinct(patient_id) %>% 
  inner_join(lung_valid_cases) %>% 
  distinct(patient_id) %>% 
  nrow()

fileConn <- file(paste0("/Shared/Statepi_Diagnosis/projects/sarcoid/","readme_sarcoid_validation.txt"))
writeLines(c("Total Number of Index Cases",
             paste0(" - Sarcoid cases with ",enroll_min," days of prior enrollment: ",n_enroll),
             paste0(" - Sarcoid cases with lung validation codes: ",n_lung),
             paste0(" - Sarcoid cases with skin validation codes: ",n_skin),
             paste0(" - Sarcoid cases with both skin and lung validation: ",n_overlap)), fileConn)
close(fileConn)


