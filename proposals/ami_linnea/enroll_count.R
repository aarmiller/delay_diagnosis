
library(tidyverse)
library(bit64)

load("/Shared/AML/truven_extracts/dx/ami/ami_index_dx_dates.RData")
index_dx_date


load("/Shared/AML/truven_extracts/dx/ami/ami_dx_visits.RData")

index <- bind_rows(all_inpatient_visits %>% 
            filter(dx_ver==10) %>% 
            distinct(enrolid,date=admdate),
          all_outpatient_visits %>% 
            filter(dx_ver==10) %>% 
            distinct(enrolid,date=svcdate),
          all_facility_visits %>% 
            filter(dx_ver==10) %>% 
            distinct(enrolid,date=svcdate),
          inpatient_services_visits %>% 
            filter(dx_ver==10) %>% 
            distinct(enrolid,date=svcdate)) %>% 
  arrange(enrolid,date) %>% 
  mutate(index = ifelse(enrolid!=lag(enrolid),1L,0L)) %>% 
  mutate(index = replace_na(index,1L)) %>% 
  filter(index==1) %>% 
  distinct(enrolid,date)


enroll_db <- src_sqlite("/Shared/Statepi_Marketscan/databases/Truven/enrollment_dbs/all_enroll.db")

ccae <- enroll_db %>% 
  tbl("ccae_collapse_enroll") %>% 
  collect() %>% 
  mutate(db = "ccae")

mdcr <- enroll_db %>% 
  tbl("mdcr_collapse_enroll") %>% 
  collect() %>% 
  mutate(db = "mdcr")

med <- enroll_db %>% 
  tbl("medicaid_collapse_enroll") %>% 
  collect() %>% 
  mutate(db = "med")

collapse_enroll <- bind_rows(ccae,mdcr,med)

rm(ccae,med,mdcr)

collapse_enroll <- inner_join(collapse_enroll,index)

tmp <- collapse_enroll %>%
  filter(date>=dtstart,date<=dtend) %>% 
  mutate(enroll_before = date-dtstart,
         enroll_after = dtend-date)


tmp %>% 
  summarise(before365 = sum(enroll_before>=365),
            before180 = sum(enroll_before>=180),
            after365 = sum(enroll_after>=365),
            after180 = sum(enroll_after>=180),
            both365 = sum(enroll_before>=365 & enroll_after>=365),
            both180 = sum(enroll_before>=180 & enroll_after>=180))



