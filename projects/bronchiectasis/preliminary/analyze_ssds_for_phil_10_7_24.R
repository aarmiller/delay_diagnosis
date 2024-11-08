
rm(list = ls())

library(tidyverse)
library(icd)
library(codeBuildr)


db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/bronchiectasis/bronchiectasis.db")

all_dx_visits <- db %>% tbl("all_dx_visits") %>% collect()

save(all_dx_visits, file = "/Shared/AML/truven_extracts/small_dbs/bronchiectasis/all_dx_visits.RData")


load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/delay_results/all_dx_visits.RData")

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/delay_results/tmp/all_dx_visits.RData")


ssd_set <- list(list(condition = "Cough",
                     dx9 = codeBuildr::children_safe("7862"),
                     dx10 = codeBuildr::children_safe("R05")),
                list(condition = "Shortness of Breath",
                     dx9 = codeBuildr::children_safe("78605"),
                     dx10 = codeBuildr::children_safe("R0602")),
                list(condition = "Pneumonia",
                     dx9 = codeBuildr::children_safe("486"),
                     dx10 = codeBuildr::children_safe("J189")))


ssd_set <- bind_rows(tibble(condition=unlist(map(ssd_set, ~.$condition)),
                 dx = unlist(map(ssd_set, ~.$dx9)),
                 dx_ver = 9L),
          tibble(condition=unlist(map(ssd_set, ~.$condition)),
                 dx = unlist(map(ssd_set, ~.$dx10)),
                 dx_ver = 10L))



prior_period <- 365*4


include_ids <- index_dx_dates %>% 
  filter(time_before_index>=prior_period) %>% 
  distinct(patient_id)

nrow(include_ids)

tmp_dx_visits <- all_dx_visits %>% 
  filter(between(days_since_index,-prior_period,-1)) %>% 
  inner_join(include_ids)

all_dx_counts <- tmp_dx_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) 

ssd_counts1 <- tmp_dx_visits %>% 
  inner_join(ssd_set) %>% 
  distinct(patient_id,days_since_index,condition) %>% 
  count(days_since_index,condition)

ssd_counts2 <- tmp_dx_visits %>% 
  inner_join(ssd_set) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  mutate(condition = "All SSDs")

all_dx_counts %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  scale_x_reverse() +
  ylab("Number of Patients") +
  xlab("Days Prior to Index") +
  theme_bw() +
  ggtitle("All Visits Prior to Diagnosis")

bind_rows(ssd_counts1,ssd_counts2) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  scale_x_reverse() +
  ylab("Number of Patients") +
  xlab("Days Prior to Index") +
  facet_wrap(~condition, scale = "free_y") +
  theme_bw() +
  ggtitle("SSD Visits Prior to Diagnosis")
