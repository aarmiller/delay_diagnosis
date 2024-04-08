

library(tidyverse)


load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/index_cases.RData")

index_cases

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/sepsis_revised10/sepsis_revised10.db")

dx_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(between(days_since_index,-180,-1)) %>% 
  filter(dx_ver==10) %>% 
  collect()