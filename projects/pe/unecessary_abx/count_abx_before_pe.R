rm(list = ls())

library(tidyverse)
library(readxl)
library(smallDB)

################
#### Params ####
################

cond_name <- "pe"

weekly_cp <- 15

upper_bound <- 365



###################
#### Load Data ####
###################

### Collect ABX visits on HPC --------------------------------------------------

abx_codes <- codeBuildr::load_rx_codes("all_abx") %>% 
  enframe() %>% 
  unnest(value) %>% 
  rename(ndc = value)


collect_codes <- distinct(abx_codes, ndcnum = ndc)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pe/pe.db")

abx_vis <- db %>% 
  tbl("all_rx_visits") %>% 
  inner_join(collect_codes, copy = TRUE) %>% 
  collect()

index_dx_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

save(abx_vis,index_dx_dates,file = "/Shared/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")

load("/Volumes/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")


### Run locally -------------------
load("/Volumes/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")
# save(abx_vis,index_dx_dates,file = "/Shared/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")


abx_counts <- rename(abx_codes,ndcnum=ndc) %>% 
  inner_join(abx_vis) %>% 
  inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(name,patient_id,days_since_index) %>% 
  count(name,days_since_index)


abx_counts %>% 
  filter(days_since_index >= -98) %>% 
  group_by(name) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  filter(name != "uncategorized") %>% 
  slice(1:6) %>% 
  select(name) %>% 
  inner_join(abx_counts) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", nrow = 3) +
  theme_bw()

all_abx_counts <- abx_counts %>% 
  filter(days_since_index >= -98) %>% 
  group_by(name) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n))


pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/pe/excess_abx/results/abx_plots/all_abx.pdf",onefile = TRUE)
for (i in 1:199){
  tmp_conds <- all_abx_counts %>%
    slice((i*6-5):(i*6)) %>% 
    select(name)
  
  p <- tmp_conds %>% 
    inner_join(abx_counts) %>% 
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    theme_bw()
  
  print(p)
}
dev.off()

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/pe/excess_abx/results/abx_plots/all_abx_week.pdf",onefile = TRUE)
for (i in 1:199){
  tmp_conds <- all_abx_counts %>%
    slice((i*6-5):(i*6)) %>% 
    select(name)
  
  p <- tmp_conds %>% 
    inner_join(abx_counts,by = join_by(name)) %>% 
    mutate(week = -(days_since_index %/% 7)) %>% 
    filter(week<53) %>% 
    group_by(name,week) %>% 
    summarise(n = sum(n)) %>% 
    ungroup() %>% 
    ggplot(aes(-week,n)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    theme_bw()
  
  print(p)
}
dev.off()


abx_counts %>% 
  mutate(week = -(days_since_index %/% 7)) %>% 
  group_by(name,week) %>% 
  summarise(n = sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(week,n)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", nrow = 3) +
  theme_bw()

abx_counts %>% 
  filter(days_since_index >= -98) %>% 
  group_by(name) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/pe/excess_abx/results/abx_before_pe.csv")


