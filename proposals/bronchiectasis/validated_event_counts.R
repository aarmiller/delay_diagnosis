

codeBuildr::load_disease_codes("bronchiectasis")


## HI Phil---possible for you to grab some prelim numbers for bronch cohort?  
## Use the below and let me know how many you find and over what time period?

## ICD-9-CM code for bronchiectasis (494.0, 494.1).  We often make it necessary 
## that at least one of these codes was given by a pulmonologist.

## Bronchiectasis diagnosis: ICD-10-CM J47.0, J47.1, J47.9  (and exclude anyone 
## with cystic fibrosis ICD-10-CM E84.0. E84.1, E84.11, E84.19, E84.8, E84.9)
## Same as above, could look also with rule mandating at least one code given by 
## a pulmonologist

## At  least two of the above ICD=9 or ICD-10 codes given within a 12 month 
## period, with at least one being from a pulmonary specialist.

rm(list = ls())

library(tidyverse)

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/bronchiectasis/bronchiectasis.db")


## Load enrolid crosswalk
load("~/Data/MarketScan/truven_extracts/small_dbs/bronchiectasis/bronchiectasis_enrolid_crosswalk.RData")
bronchiectasis_crosswalk <- enrolids

tmp <- codeBuildr::load_disease_codes("bronchiectasis")$bronchiectasis

bronchiectasis_codes <- bind_rows(tibble(dx=tmp$icd9_codes,
                                         dx_ver = 9L),
                                  tibble(dx=tmp$icd10_codes,
                                         dx_ver = 10L))
          

bronchiectasis_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(bronchiectasis_codes, copy = TRUE) %>% 
  collect()

collect_pulmonologist_vis <- function(year){
  
  tmp_ccae <- db %>% 
    tbl(paste0("outpatient_core_ccae_",year)) %>% 
    filter(stdprov %in% c("295","425")) %>% 
    distinct(patient_id,svcdate) %>% 
    collect()
  
  tmp_mdcr <- db %>% 
    tbl(paste0("outpatient_core_mdcr_",year)) %>% 
    filter(stdprov %in% c("295","425")) %>% 
    distinct(patient_id,svcdate) %>% 
    collect()
  
  res <- bind_rows(tmp_ccae,tmp_mdcr)
  
  if (year %in% 14:21){
    
    tmp_med <- db %>% 
      tbl(paste0("outpatient_core_medicaid_",year)) %>% 
      filter(stdprov %in% c("295","425")) %>% 
      distinct(patient_id,svcdate) %>% 
      collect()
    
    return(bind_rows(res,tmp_med))
    
  } else {
    
    return(res)
    
  }
  
}


pulm_vis_out <- tibble(year = str_pad(1:22, pad = 0, width = 2)) %>% 
  mutate(data = map(year,collect_pulmonologist_vis)) %>% 
  unnest(data) 

pulm_vis_out <- pulm_vis_out %>% 
  distinct(patient_id,svcdate)

bronchiectasis_visits %>% 
  distinct(patient_id)

out_pulm_dx_cases <- bronchiectasis_visits %>% 
  distinct(patient_id,svcdate=date) %>% 
  inner_join(pulm_vis_out) %>% 
  distinct(patient_id)


inpatient_pulm <- db %>% 
  tbl("inpatient_services_core") %>% 
  filter(stdprov %in% c("295","425")) %>% 
  distinct(patient_id,caseid,admdate,disdate,svcdate) %>% 
  collect()

hosp_pulm_dx_cases <- bronchiectasis_visits %>% 
  filter(inpatient==1) %>% 
  distinct(patient_id,date) %>% 
  inner_join(inpatient_pulm) %>% 
  filter(date == svcdate | (date<=disdate & date>=admdate)) %>% 
  distinct(patient_id)

# FINAL PULM ID CASES
pulm_id_cases <- bind_rows(out_pulm_dx_cases,hosp_pulm_dx_cases) %>% 
  distinct(patient_id)



### Remove CF -------------------

load("~/Data/MarketScan/truven_extracts/small_dbs/cf/cf_enrolid_crosswalk.RData")
cf_enrolids <- enrolids %>% distinct(enrolid)


cf_ids <- bronchiectasis_crosswalk %>% 
  inner_join(cf_enrolids) %>% 
  distinct(patient_id)


final_cases <- pulm_id_cases %>% 
  anti_join(cf_ids)


## Cases assigned by specialist and 2 codes within 12mo window

twodx_12mo <- bronchiectasis_visits %>% 
  distinct(patient_id,date) %>% 
  arrange(patient_id,date) %>% 
  mutate(days_since = date-lag(date)) %>% 
  mutate(days_since = ifelse(patient_id!=lag(patient_id), NA,days_since)) %>% 
  filter(!is.na(days_since)) %>% 
  filter(days_since<=365) %>% 
  distinct(patient_id)

final_cases %>% 
  inner_join(twodx_12mo)


### Compute enrollment time prior --------

index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

enroll_scale <- 4

index_dx_dates %>% 
  filter(time_before_index>=365*enroll_scale) %>% 
  nrow()

index_dx_dates %>% 
  anti_join(cf_ids,by = join_by(patient_id)) %>% 
  filter(time_before_index>=365*enroll_scale) %>% 
  nrow()

index_dx_dates %>% 
  anti_join(cf_ids,by = join_by(patient_id)) %>% 
  inner_join(pulm_id_cases,by = join_by(patient_id)) %>% 
  filter(time_before_index>=365*enroll_scale) %>% 
  nrow()

index_dx_dates %>% 
  anti_join(cf_ids,by = join_by(patient_id)) %>% 
  inner_join(pulm_id_cases,by = join_by(patient_id)) %>% 
  inner_join(twodx_12mo,by = join_by(patient_id)) %>% 
  filter(time_before_index>=365*enroll_scale) %>% 
  nrow()

