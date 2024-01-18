
library(tidyverse)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis.db")

# Load Index Visits and Crosswalk
enrolid_crosswalk <- db %>% tbl("enrolid_crosswalk") %>% collect()

index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

# Add family ID
enrolid_crosswalk <- enrolid_crosswalk %>% 
  mutate(efamid = enrolid %/% 100)

# Load validated ids
validated_ids <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("87798","86615")) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect() %>% 
  distinct(patient_id)

# Find validated family members
all_validated_ids <- enrolid_crosswalk %>% 
  inner_join(validated_ids) %>% 
  distinct(efamid) %>% 
  inner_join(enrolid_crosswalk)

##################
#### Lab Data ####


# labs_db <- DBI::dbConnect(RSQLite::SQLite(),"/Shared/Statepi_Marketscan/databases/Truven/lab_dbs/labs.db")
# 
# pertus_codes <- labs_db %>% 
#   tbl("loinc_labels") %>% 
#   collect() %>% 
#   filter(str_detect(long_common_name,"pertus")) %>% 
#   distinct(loinccd,long_common_name)
# 
# pull_codes <- pertus_codes$loinccd
# 
# loinc_visits <- tibble()
# 
# for (i in 16:21){
#   print(i)
#   
#   tmp <- labs_db %>% tbl(paste0("labs_mdcr_",i)) %>% filter(loinccd %in% pull_codes) %>% collect()
#   
#   loinc_visits <- bind_rows(loinc_visits,tmp) %>% distinct()
#   
#   tmp <- labs_db %>% tbl(paste0("labs_ccae_",i)) %>% filter(loinccd %in% pull_codes) %>% collect()
#   
#   loinc_visits <- bind_rows(loinc_visits,tmp) %>% distinct()
#   
# }
# 
# 
# pertussis_labs <- loinc_visits
# 
# save(pertussis_labs, file = "/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis_labs.RData")

load("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis_labs.RData")

distinct(all_validated_ids,enrolid) %>% 
  inner_join(pertussis_labs) %>% 
  distinct(enrolid)

distinct(all_validated_ids,enrolid) %>% 
  inner_join(pertussis_labs) %>% 
  filter(result>refhigh | result < reflow) %>% 
  distinct(enrolid)


loinc_visits %>% 
  count(resltcat)

loinc_visits %>% count(abnormal)


loinc_visits %>% filter(abnormal %in% c("A","H")) %>% count(resltcat)

pertus_lab_test <- loinc_visits %>% distinct(enrolid)


loinc_visits %>% 
  filter(resltcat %in% c("DET","POS"))

loinc_labels %>% filter(str_detect(long_common_name,"pertus"))

load("/Volumes/Statepi_Marketscan/databases/Truven/lab_dbs/distinct_loinc.RData")

tmp <- write_db %>% tbl("labs_mdcr_16") %>% distinct(loinccd) %>% collect()


write_db %>% tbl("labs_ccae_19") %>% filter(loinccd %in% c("23826-1","92128-8",
                                                           "11585-7","43913-3")) %>% count(resltcat)
