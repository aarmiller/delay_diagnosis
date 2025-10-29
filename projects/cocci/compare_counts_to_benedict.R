rm(list = ls())
library(tidyverse)
# db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/cocci/cocci.db")
db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/cocci/cocci.db")

tmp <- codeBuildr::load_disease_codes("cocci",return_tibble = TRUE)

library(codeBuildr)


benedict_set <- c("R634","R07","R6883","R05","R06","L52","R538","R5081","R509",
  "R61","M791","M255","J960","J985","R59","R0902","R918","J90","R911",
  "J06","J20","J09","J10","J11","C34","J984","J989","J12","J13","J14","J15",
  "J16","J17","J18","J22")

benedict_set <- children_safe(benedict_set)



index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

include_cases <- index_dates %>% 
  filter(index_date>=as.integer(ymd("2017-01-01"))) %>% 
  filter(time_before_index>=365)

medicaid_ids <- db %>% 
  tbl("medicaid_collapse_enroll") %>% 
  distinct(patient_id) %>% 
  collect()

include_cases <- include_cases %>% 
  anti_join(medicaid_ids)

benedict_ssds <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(dx_ver==10 & dx %in% benedict_set) %>% 
  filter(days_since_index>=-92 & days_since_index<=-1) %>% 
  collect() 


(benedict_ssds %>% 
  inner_join(include_cases) %>% 
  distinct(patient_id) %>% 
  nrow()) / nrow(include_cases)


our_set <- codeBuildr::load_ssd_codes("cocci") %>% 
  filter(type == "icd10") %>% 
  .$code

our_ssds <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(dx_ver==10 & dx %in% our_set) %>% 
  filter(days_since_index>=-92 & days_since_index<=-1) %>% 
  collect() 




(our_ssds %>% 
    inner_join(include_cases) %>% 
    distinct(patient_id) %>% 
    nrow()) / nrow(include_cases)


     