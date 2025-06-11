library(tidyverse)
library(lubridate)
library(haven)
library(bit64)


# Specify cond and output path -------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
cond <- args[1]
# cond <- "sepsis_kaiser"

out_path <- paste0("/Shared/AML/kaiser_data/", str_remove(cond, "_kaiser"), "/delay_data/")
ifelse(!dir.exists(out_path), dir.create(out_path), "directory exists")

# Build database ---------------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), paste0(out_path, cond, ".db"))
DBI::dbListTables(con)

# Create index_dx_dates table ---------------------------------------------------

index_dx_dates <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", str_remove(cond, "_kaiser"), "/", 
                                        str_remove(cond, "_kaiser"), "_cohort_12sep23final.sas7bdat"))

index_dx_dates <- index_dx_dates %>% 
  select(patient_id = STUDYID, index_date, time_before_index = enroll_before) %>% 
  mutate(index_date = as.integer(index_date)) %>% 
  distinct()

dplyr::copy_to(con, index_dx_dates, "index_dx_dates", temporary = FALSE, overwrite = TRUE)

# Add all_dx_visits to db ------------------------------------------------------

all_enconters <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", str_remove(cond, "_kaiser"),
                                        "/", str_remove(cond, "_kaiser"), "_enc_diag_12sep23final.sas7bdat"))

all_dx_visits <- all_enconters %>% 
  select(patient_id = STUDYID, dx, dx_ver = dx_codetype,
         encounter_date, encounter_type, encounter_subtype) %>% 
  mutate(encounter_date =as.integer(encounter_date)) %>% 
  inner_join(index_dx_dates %>% select(-time_before_index), by = c("patient_id")) %>% 
  mutate(days_since_index = encounter_date-index_date) %>% 
  select(-index_date) %>% 
  mutate(dx = str_remove_all(dx , "\\."))

dplyr::copy_to(con, all_dx_visits, "all_dx_visits", temporary = FALSE, overwrite = TRUE)


# Add all_rx_visits to db ------------------------------------------------------
med_admin <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", str_remove(cond, "_kaiser"),
                                    "/", str_remove(cond, "_kaiser"), "_med_admin_12sep23final.sas7bdat"))

all_rx_visits <- med_admin %>% 
  rename(patient_id = STUDYID,
         date = admin_date,
         ndcnum = ndc) %>% 
  mutate(date =as.integer(date))

dplyr::copy_to(con, all_rx_visits, "all_rx_visits", temporary = FALSE, overwrite = TRUE)

# Add tm_full to db -----------------------------------------------------------------

#disdate not NA for HOSPITAL AND OVERNIGHT ENCOUNTERS ONLY according to data dictionary

tm_full <- all_enconters %>% 
  select(patient_id = STUDYID, 
         case_id = encounter_id,
         admdate= encounter_date, disdate = discharge_date,
         encounter_type) %>% 
  mutate(ind = 1L) %>% 
  distinct() %>% 
  pivot_wider(id_cols = c("patient_id", "case_id", "admdate", "disdate"),
              names_from = encounter_type,
              values_from = ind) %>% 
  mutate(across(!patient_id & !case_id & !admdate & !disdate, ~ifelse(is.na(.x), 0L, .x))) %>% 
  mutate(svcdate = admdate)

tm_full <- tm_full %>% 
  mutate(inpatient = as.integer( (IP + IS) >= 1)) %>% 
  mutate(outpatient = ifelse(AV == 1, 1L, 0L)) %>% 
  mutate(other = as.integer( (LO + OE + RO + VC) >= 1)) %>% 
  rename(ed = ED) 

dplyr::copy_to(con, tm_full, "tm_full", temporary = FALSE, overwrite = TRUE)

# Add tm to db -----------------------------------------------------------------

tm <- all_dx_visits %>% 
  select(patient_id, svcdate= encounter_date, encounter_type) %>% 
  mutate(ind = 1L) %>% 
  distinct() %>% 
  pivot_wider(id_cols = c("patient_id", "svcdate"),
              names_from = encounter_type,
              values_from = ind) %>% 
  mutate(across(!patient_id & !svcdate, ~ifelse(is.na(.x), 0L, .x))) 

# setting_inds <- all_dx_visits %>% 
#   select(patient_id, svcdate= encounter_date, encounter_type, encounter_subtype) %>% 
#   mutate(setting = ifelse(encounter_type == "IP" | encounter_type == "IS", "inpatient",
#                           ifelse(encounter_type == "AV" & encounter_subtype == "OB", "obs_stay",
#                                  ifelse(encounter_type == "AV" & encounter_subtype != "OB", "outpatient" , "other")))) %>% 
#   filter(setting != "other") %>% 
#   mutate(ind = 1L) %>% 
#   distinct(patient_id, svcdate, setting, ind) %>% 
#   pivot_wider(id_cols = c("patient_id", "svcdate"),
#               names_from = setting,
#               values_from = ind) %>% 
#   mutate(across(!patient_id & !svcdate, ~ifelse(is.na(.x), 0L, .x))) 

tm <- tm %>% 
  mutate(inpatient = as.integer( (IP + IS) >= 1)) %>% 
  mutate(outpatient = ifelse(AV == 1, 1L, 0L)) %>% 
  mutate(other = as.integer( (LO + OE + RO + VC) >= 1)) %>% 
  rename(ed = ED) 
  
dplyr::copy_to(con, tm, "tm", temporary = FALSE, overwrite = TRUE)

# Add dx_dates_alg1 to db ------------------------------------------------------

dx_dates_alg1 <- all_dx_visits %>% 
  mutate(inpatient = ifelse(encounter_type == "IP" | encounter_type == "IS", 1L, 0L)) %>% 
  select(patient_id, date = encounter_date, dx, dx_ver, inpatient) %>% 
  distinct()

dplyr::copy_to(con, dx_dates_alg1, "dx_dates_alg1", temporary = FALSE, overwrite = TRUE)

# Add stdplac_visits to db ------------------------------------------------------

# essentially just more granular encounter type (i.e., ENCOUNTER_SUBTYPE)

stdplac_visits <- all_dx_visits %>% 
  select(patient_id, svcdate = encounter_date, encounter_subtype) %>% 
  distinct()

dplyr::copy_to(con, stdplac_visits, "stdplac_visits", temporary = FALSE, overwrite = TRUE)

# Add stdprov_visits to db -----------------------------------------------------

# Cannot do this as we dont have provider information

# Add svcscat_visits to db -----------------------------------------------------

# Cannot do this as we dont have this information

# Close connection -------------------------------------------------------------

DBI::dbListTables(con)
DBI::dbDisconnect(con)
