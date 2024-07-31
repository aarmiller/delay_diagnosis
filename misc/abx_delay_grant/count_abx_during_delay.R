

library(tidyverse)

rm(list = ls())

cond_name <- "ami"

load("/Shared/AML/params/delay_any_params.RData")
delay_params <- delay_any_params[[cond_name]]
rm(delay_any_params)

delay_base_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")

con <- DBI::dbConnect(RSQLite::SQLite(), paste0("/Shared/AML/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

### Pull Index Visits ###

index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=delay_params$upper_bound)


### Pull Antibiotics ###

abx_codes <- unlist(codeBuildr::load_rx_codes("all_abx"),use.names = F)

abx_vis <- con %>% 
  tbl("all_rx_visits") %>% 
  collect() %>% 
  filter(ndcnum %in% abx_codes)


#### Count patients with antibiotic before diagnosis ####

index_count <- distinct(index_dx_dates,patient_id) %>% nrow()

abx_count <- index_dx_dates %>% 
  select(patient_id,index_date) %>% 
  inner_join(abx_vis) %>% 
  filter(date<index_date & date>=(index_date-delay_params$cp)) %>% 
  distinct(patient_id) %>% 
  nrow()

round(100*abx_count/index_count,2)

##########################################################
### Count antibiotics during diagnostic opportunities ####
##########################################################

load(paste0(delay_base_path,"ssd_visit/sim_res.RData"))

abx_inds <- abx_vis %>% 
  inner_join(index_dx_dates) %>% 
  mutate(days_since_index = date-index_date) %>% 
  distinct(patient_id,days_since_index) %>% 
  filter(between(days_since_index,-delay_params$cp,-1)) %>% 
  mutate(abx_vis = 1L)

### SSD Visits - Individual Opportunity window ---------------------------------


tmp <- sim_res %>% 
  inner_join(sim_res_sim_obs)

tmp <- tmp %>% 
  group_by(trial,patient_id) %>% 
  summarise(delay = -min(days_since_index))

abx_counts <- tmp %>% 
  ungroup() %>% 
  inner_join(abx_inds) %>% 
  filter(-days_since_index<delay) %>% 
  distinct(trial,patient_id) %>% 
  count(trial, name = "abx_count")

miss_count1 <- tmp %>% 
  ungroup() %>% 
  count(trial, name = "miss_count") %>% 
  left_join(abx_counts) %>% 
  mutate(delay_pct = 100*abx_count/miss_count) %>% 
  summarise(mean_pct = mean(delay_pct),
            lo_pct = quantile(delay_pct, probs = c(0.025)),
            hi_pct = quantile(delay_pct, probs = c(0.975)))

### SSD Visits - Aggregate Opportunity window ----------------------------------


abx_counts <- tmp %>% 
  ungroup() %>% 
  inner_join(distinct(abx_inds,patient_id)) %>% 
  count(trial, name = "abx_count")

miss_count2 <- tmp %>% 
  ungroup() %>% 
  count(trial, name = "miss_count") %>% 
  left_join(abx_counts) %>% 
  mutate(delay_pct = 100*abx_count/miss_count) %>% 
  summarise(mean_pct = mean(delay_pct),
            lo_pct = quantile(delay_pct, probs = c(0.025)),
            hi_pct = quantile(delay_pct, probs = c(0.975)))


### All Visits - Individual Opportunity window ---------------------------------

load(paste0(delay_base_path,"any_visit/sim_res.RData"))

tmp <- sim_res %>% 
  inner_join(sim_res_sim_obs)

tmp <- tmp %>% 
  group_by(trial,patient_id) %>% 
  summarise(delay = -min(days_since_index))

abx_counts <- tmp %>% 
  ungroup() %>% 
  inner_join(abx_inds) %>% 
  filter(-days_since_index<delay) %>% 
  distinct(trial,patient_id) %>% 
  count(trial, name = "abx_count")

miss_count3 <- tmp %>% 
  ungroup() %>% 
  count(trial, name = "miss_count") %>% 
  left_join(abx_counts) %>% 
  mutate(delay_pct = 100*abx_count/miss_count) %>% 
  summarise(mean_pct = mean(delay_pct),
            lo_pct = quantile(delay_pct, probs = c(0.025)),
            hi_pct = quantile(delay_pct, probs = c(0.975)))


### All Visits - Aggregate Opportunity window ----------------------------------


abx_counts <- tmp %>% 
  ungroup() %>% 
  inner_join(distinct(abx_inds,patient_id)) %>% 
  count(trial, name = "abx_count")

miss_count4 <- tmp %>% 
  ungroup() %>% 
  count(trial, name = "miss_count") %>% 
  left_join(abx_counts) %>% 
  mutate(delay_pct = 100*abx_count/miss_count) %>% 
  summarise(mean_pct = mean(delay_pct),
            lo_pct = quantile(delay_pct, probs = c(0.025)),
            hi_pct = quantile(delay_pct, probs = c(0.975)))


###########################
#### Aggregate Results ####
###########################

abx_miss_count <- bind_rows(mutate(miss_count1,group = "SSD - Individual"),
                            mutate(miss_count2,group = "SSD - Aggregate"),
                            mutate(miss_count3,group = "ANY - Individual"),
                            mutate(miss_count4,group = "ANY - Aggregate")) %>% 
  mutate_at(vars(mean_pct:hi_pct),~round(.,1)) %>% 
  mutate(abx_rate = paste0(mean_pct," (",lo_pct,"-",hi_pct,")")) %>% 
  select(group,abx_rate)

abx_miss_count %>% 
  write_csv(paste0("/Shared/AML/grant_proposal_work/abx_delay_2024/data/abx_delay_counts/",cond_name,".csv"))
