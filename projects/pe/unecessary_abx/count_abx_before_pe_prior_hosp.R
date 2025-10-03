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

# find inpatient visit dates
in_tm1 <- db %>% 
  tbl("tm") %>% 
  filter(inpatient==1) %>% 
  select(patient_id,svcdate) %>% 
  collect()

in_tm2 <- db %>% 
  tbl("tm_full") %>% 
  filter(inpatient==1) %>% 
  select(patient_id,svcdate) %>% 
  collect()

inpatient_dates <- in_tm1 %>% 
  full_join(in_tm2) %>% 
  distinct()


save(abx_vis,index_dx_dates,inpatient_dates,file = "/Shared/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")


### Run locally -------------------
# load("/Volumes/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")
# save(abx_vis,index_dx_dates,file = "/Shared/AML/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")

rm(list = ls())
load("~/Data/MarketScan/truven_extracts/small_dbs/pe/tmp_transfer/abx_visits.RData")

## Params ##
inpatient_window <- 90

############

abx_codes <- codeBuildr::load_rx_codes("all_abx") %>% 
  enframe() %>% 
  unnest(value) %>% 
  rename(ndc = value)


prior_in_inds <- inpatient_dates %>% 
  inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
  mutate(days_since_index = svcdate-index_date) %>% 
  filter(between(days_since_index,-inpatient_window,-1)) %>% 
  distinct(patient_id) %>% 
  mutate(prior_inpatient = 1L)

abx_counts <- rename(abx_codes,ndcnum=ndc) %>% 
  inner_join(abx_vis) %>% 
  inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  left_join(prior_in_inds) %>% 
  mutate(prior_inpatient = replace_na(prior_inpatient,0L)) %>% 
  distinct(name,patient_id,days_since_index,prior_inpatient) %>% 
  count(name,prior_inpatient,days_since_index)


library("ggh4x")
install.packages("ggh4x")

abx_counts %>% 
  filter(days_since_index >= -98) %>% 
  group_by(name) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  filter(name != "uncategorized") %>% 
  slice(1:3) %>% 
  select(name) %>% 
  inner_join(abx_counts) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_grid(name~prior_inpatient, scales = "free") +
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


