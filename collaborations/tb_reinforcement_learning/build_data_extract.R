
library(tidyverse)

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/tb/tb.db")


## Timemap ---------------------------------------------------------------------

timemap <- db %>% 
  tbl("tm") %>% 
  select(patient_id,date=svcdate,mdcr,ccae,medicaid,outpatient,ed,obs_stay,inpatient,rx) %>% 
  collect()

timemap %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/timemap.csv")

rm(timemap)

## Index DX --------------------------------------------------------------------

index_dx <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dx %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/index_tb_date.csv")

rm(index_dx)

## TB DX -----------------------------------------------------------------------

codeBuildr::load_disease_codes("tb",return_tibble = TRUE)

tb_dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(codeBuildr::load_disease_codes("tb",return_tibble = TRUE), 
             copy = TRUE,
             by = join_by(dx, dx_ver)) %>% 
  collect()

tb_dx_visits <- tb_dx_visits %>% 
  distinct(patient_id,dx,dx_ver,inpatient,date)

tb_dx_visits %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/tb_dx_visits.csv")

rm(tb_dx_visits)
  

## Enrollment ------------------------------------------------------------------

# Collapsed enrollment
tmp1 <- db %>% 
  tbl("ccae_mdcr_collapse_enroll") %>% 
  collect()

tmp2 <- db %>% 
  tbl("medicaid_collapse_enroll") %>% 
  collect()

collapse_enroll <- bind_rows(tmp1,tmp2)

collapse_enroll %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/enrollment_periods.csv")

# Enrollment detail
rm(collapse_enroll,tmp1,tmp2)

all_enroll_ccae_mdcr <- db %>% 
  tbl("all_enroll") %>% 
  collect()

all_enroll_medicaid <- db %>% 
  tbl("all_enroll_medicaid") %>% 
  collect()

all_enroll_ccae_mdcr %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/all_enroll_ccae_mdcr.csv")

all_enroll_medicaid %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/all_enroll_medicaid.csv")

rm(all_enroll_ccae_mdcr,all_enroll_medicaid)


## Diagnoses -------------------------------------------------------------------

all_dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  select(patient_id:date) %>% 
  collect()

all_dx_visits %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/all_dx_visits.csv")

rm(all_dx_visits)

codeBuildr::all_icd_labels %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/icd_labels.csv")

## Procedures ------------------------------------------------------------------

all_proc_visits <- db %>% 
  tbl("all_proc_visits") %>% 
  select(patient_id:date) %>% 
  collect()

all_proc_visits %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/all_proc_visits.csv")

rm(all_proc_visits)

## Medications -----------------------------------------------------------------

all_rx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  collect()

all_rx_visits %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/all_rx_visits.csv")

redbook <- haven::read_sas("~/Data/MarketScan/redbook/redbook_Jan_2025.sas7bdat")

names(redbook) <- tolower(names(redbook))

redbook %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/redbook.csv")

rm(all_rx_visits,redbook)

## Demographics ----------------------------------------------------------------

## Providers -------------------------------------------------------------------

stdprov_dates <- db %>% 
  tbl("stdprov_visits") %>% 
  select(patient_id,date=svcdate,stdprov) %>% 
  collect()

stdprov_dates %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/stdprov_dates.csv")

rm(stdprov_dates)

## Places of Care --------------------------------------------------------------

stdplac_dates <- db %>% 
  tbl("stdplac_visits") %>% 
  select(patient_id,date=svcdate,stdplac) %>% 
  collect()

stdplac_dates %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/stdplac_dates.csv")

rm(stdplac_dates)

## States ----------------------------------------------------------------------

collect_location <- function(year){
  tmp1 <- db %>% 
    tbl(paste0("enrollment_detail_ccae_",year)) %>% 
    select(patient_id,month,dtend,dtstart,msa,egeoloc) %>% 
    mutate(msa = as.character(msa)) %>% 
    collect()
  
  tmp2 <- db %>% 
    tbl(paste0("enrollment_detail_mdcr_",year)) %>% 
    select(patient_id,month,dtend,dtstart,msa,egeoloc) %>% 
    mutate(msa = as.character(msa)) %>% 
    collect()
  
  bind_rows(tmp1,tmp2) %>% 
    mutate(year = as.integer(year)+2000L) %>% 
    select(patient_id,year,month:egeoloc)
}

tmp <- tibble(year = str_pad(1:23,width = 2,pad = 0)) %>% 
  mutate(data = map(year,collect_location))

patient_location <- tmp %>% 
  select(data) %>% 
  unnest(data)

patient_location %>% 
  write_csv("/Volumes/Statepi_Shared/users/aarmille/tb_delay_learning/data/patient_location.csv")

rm(list = ls())

