
library(tidyverse)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/sarcoid/sarcoid.db")


# Procedure codes for identifying Sarcoid Lung
lung_procs <- c("31628","31624","31629","3324","31622","31620","31623","31633",
                 "31632","3323","3322","31645","31625","88172","3324","10022",
                 "76942","32405","38510","38505","3328","38525","3326","4011")

# Procedure codes for identifying Sarcoid Skin
skin_procs <- c("11100", "11101")

# Pull procedure codes
procs <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c(skin_procs,lung_procs)) %>% 
  filter(between(days_since_index,-30,1)) %>% 
  collect()


sarcoid_skin <- procs %>% 
  filter(proc %in% skin_procs) %>% 
  distinct(patient_id)


sarcoid_lung <- procs %>% 
  filter(proc %in% lung_procs) %>% 
  distinct(patient_id)


save(sarcoid_lung, file = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/sarcoid_lung_patids.RData")
save(sarcoid_skin, file = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_skin/sarcoid_skin_patids.RData")
