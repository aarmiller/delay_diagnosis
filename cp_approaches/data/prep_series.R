
library(tidyverse)

rm(list = ls())

load("/Volumes/AML/params/final_delay_params.RData")

cond_name <- "sarcoid_skin"

delay_params <- final_delay_params[[cond_name]]

ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>% 
  filter(type %in% c("icd9","icd10")) %>% 
  mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
  select(dx = code,dx_ver)

load(paste0("/Volumes/Statepi_Diagnosis/projects/sarcoid/",cond_name,"/index_cases.RData"))
load(paste0("/Volumes/Statepi_Diagnosis/prelim_results/sarcoid/delay_results/all_dx_visits.RData"))

index_cases

patient_ids <- index_cases %>% distinct(patient_id)
all_dx_visits <- all_dx_visits %>% inner_join(patient_ids)

visit_counts <- all_dx_visits %>%
  distinct(patient_id,dx_ver,days_since_index) %>%
  count(dx_ver,days_since_index)

# populate missing valules in visit counts (i.e., assign 0 to days missing)
visit_counts <- tibble(days_since_index=-delay_params$upper_bound:delay_params$upper_bound) %>%
  mutate(dx_ver=map(days_since_index,~c(9,10))) %>%
  unnest(dx_ver) %>%
  arrange(dx_ver,days_since_index) %>%
  left_join(visit_counts,by = c("days_since_index", "dx_ver")) %>%
  mutate(n = replace_na(n,0))

tmp <- all_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index)

visit_counts <- bind_rows(tmp,visit_counts %>%
                            filter(days_since_index<=0))

count_data_ssd <- all_dx_visits %>%
  inner_join(ssd_codes, by = "dx") %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  left_join(tibble(days_since_index=min(visit_counts$days_since_index):-1),., by = "days_since_index") %>% # in case there are 0 days
  mutate(n = replace_na(n,0L)) %>% 
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))


count_data_ssd %>% 
  ggplot(aes(days_since_index,n)) +
  geom_point()

count_data_ssd <- count_data_ssd %>% mutate(dow = paste0("dow_",dow))


write_csv(count_data_ssd,file = paste0("cp_approaches/data/",cond_name,"_ssd_counts.csv"))

