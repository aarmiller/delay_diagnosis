

library(tidyverse)
library(readxl)

upper_bound <- 365

################
#### Blasto ####
################

cond_name <- "blasto"

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx") %>% 
  rename(name = `Antibiotic Name`) %>% 
  gather(key = disease, value = include, -name)

load("/Volumes/argon_home/projects/delay_diagnosis/excess_abx/data/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == "sporo") %>% 
  filter(include %in% c("Y","?")) %>% 
  filter(!(name %in% c("linezolid","cefdinir","dicloxacillin"))) %>%
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

#### Load Disease Data ---------------------------------------------------------

db <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

abx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  filter(ndcnum %in% ndc_codes) %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dates <- index_dates %>% 
  filter(time_before_index>=upper_bound) %>% 
  distinct(patient_id,index_date)

# filter out 2023
index_dates <- index_dates %>% 
  filter(index_date<as.integer(ymd("2023-01-01"))) 

## Pull rx enroll locally for ccae/mdcr ----------------------------------------

tmp_collect_enroll <- function(source,year){
  db %>% 
    tbl(paste0("enrollment_detail_",source,"_",year)) %>% 
    filter(rx==1) %>% 
    select(patient_id,dtstart,dtend) %>% 
    collect() 
}

tmp_enroll_info <- collect_plan(db) %>% 
  filter(source!="medicaid") %>% 
  mutate(data = map2(source,year,tmp_collect_enroll))

load("/Volumes/AML/tmp_transfer/excess_abx/blasto_med_rx_enroll.RData")

# assemble rx enrollment dates
rx_enroll <- bind_rows(tmp_enroll_info %>% 
                         unnest(data) %>% 
                         select(patient_id,dtstart,dtend),
                       select(med_rx_enroll,patient_id,dtstart,dtend)) %>% 
  arrange(patient_id,dtstart) %>% 
  mutate(gap = dtstart != lag(dtend)+1) %>% 
  mutate(gap = replace_na(gap,FALSE)) %>% 
  mutate(new_id = patient_id!=lag(patient_id)) %>% 
  mutate(new_id = replace_na(new_id,FALSE)) %>% 
  mutate(tmp =  ifelse(gap==TRUE | new_id == TRUE, 1L, 0L)) %>% 
  mutate(period = cumsum(tmp)) %>% 
  group_by(patient_id,period) %>% 
  summarise(dtstart = min(dtstart),
            dtend = max(dtend)) %>% 
  ungroup()

# Final index dates
index_dates <- index_dates %>% 
  inner_join(rx_enroll) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  mutate(rx_days_before = index_date-dtstart) %>% 
  filter(rx_days_before>=365) %>% 
  select(patient_id,index_date) 


#### Assemble Final Data -------------------------------------------------------

# distinct abx visit days (to use for analysis)
final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,upper_bound)) %>% 
  group_by(patient_id,name,days_since_index) %>% 
  summarise(daysupp = sum(daysupp)) %>% 
  ungroup() %>% 
  mutate(period = days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

final_abx_counts <- final_abx_visits %>% 
  count(name,week)

# remove incomplete weeks
week_exclude <- final_abx_visits %>% 
  distinct(week,period) %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)

final_abx_counts <- final_abx_counts %>% 
  anti_join(week_exclude)

final_abx_visits <- final_abx_visits %>% 
  anti_join(week_exclude)

out_data <- final_abx_visits %>% 
  distinct(patient_id,name,week) %>% 
  mutate(disease = "Blastomycosis")

rm(list = ls()[!(ls() %in% c("out_data","upper_bound"))])

###############
#### Cocci ####
###############

cond_name <- "cocci"

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx") %>% 
  rename(name = `Antibiotic Name`) %>% 
  gather(key = disease, value = include, -name)

load("/Volumes/argon_home/projects/delay_diagnosis/excess_abx/data/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == "sporo") %>% 
  filter(include %in% c("Y","?")) %>% 
  filter(!(name %in% c("linezolid","cefdinir","dicloxacillin"))) %>%
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

#### Load Disease Data ---------------------------------------------------------

db <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

abx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  filter(ndcnum %in% ndc_codes) %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dates <- index_dates %>% 
  filter(time_before_index>=upper_bound) %>% 
  distinct(patient_id,index_date)

# filter out 2023
index_dates <- index_dates %>% 
  filter(index_date<as.integer(ymd("2023-01-01"))) 

## Pull rx enroll locally for ccae/mdcr ----------------------------------------

tmp_collect_enroll <- function(source,year){
  db %>% 
    tbl(paste0("enrollment_detail_",source,"_",year)) %>% 
    filter(rx==1) %>% 
    select(patient_id,dtstart,dtend) %>% 
    collect() 
}

tmp_enroll_info <- collect_plan(db) %>% 
  filter(source!="medicaid") %>% 
  mutate(data = map2(source,year,tmp_collect_enroll))

load("/Volumes/AML/tmp_transfer/excess_abx/cocci_med_rx_enroll.RData")

# assemble rx enrollment dates
rx_enroll <- bind_rows(tmp_enroll_info %>% 
                         unnest(data) %>% 
                         select(patient_id,dtstart,dtend),
                       select(med_rx_enroll,patient_id,dtstart,dtend)) %>% 
  arrange(patient_id,dtstart) %>% 
  mutate(gap = dtstart != lag(dtend)+1) %>% 
  mutate(gap = replace_na(gap,FALSE)) %>% 
  mutate(new_id = patient_id!=lag(patient_id)) %>% 
  mutate(new_id = replace_na(new_id,FALSE)) %>% 
  mutate(tmp =  ifelse(gap==TRUE | new_id == TRUE, 1L, 0L)) %>% 
  mutate(period = cumsum(tmp)) %>% 
  group_by(patient_id,period) %>% 
  summarise(dtstart = min(dtstart),
            dtend = max(dtend)) %>% 
  ungroup()

# Final index dates
index_dates <- index_dates %>% 
  inner_join(rx_enroll) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  mutate(rx_days_before = index_date-dtstart) %>% 
  filter(rx_days_before>=365) %>% 
  select(patient_id,index_date) 


#### Assemble Final Data -------------------------------------------------------

# distinct abx visit days (to use for analysis)
final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,upper_bound)) %>% 
  group_by(patient_id,name,days_since_index) %>% 
  summarise(daysupp = sum(daysupp)) %>% 
  ungroup() %>% 
  mutate(period = days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

final_abx_counts <- final_abx_visits %>% 
  count(name,week)

# remove incomplete weeks
week_exclude <- final_abx_visits %>% 
  distinct(week,period) %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)

final_abx_counts <- final_abx_counts %>% 
  anti_join(week_exclude)

final_abx_visits <- final_abx_visits %>% 
  anti_join(week_exclude)

tmp_out <- final_abx_visits %>% 
  distinct(patient_id,name,week) %>% 
  mutate(disease = "Coccidioidomycosis")

out_data <- bind_rows(out_data,tmp_out)

rm(list = ls()[!(ls() %in% c("out_data","upper_bound"))])

###############
#### Histo ####
###############

cond_name <- "histo"

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx") %>% 
  rename(name = `Antibiotic Name`) %>% 
  gather(key = disease, value = include, -name)

load("/Volumes/argon_home/projects/delay_diagnosis/excess_abx/data/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == "sporo") %>% 
  filter(include %in% c("Y","?")) %>% 
  filter(!(name %in% c("linezolid","cefdinir","dicloxacillin"))) %>%
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

#### Load Disease Data ---------------------------------------------------------

db <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

abx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  filter(ndcnum %in% ndc_codes) %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dates <- index_dates %>% 
  filter(time_before_index>=upper_bound) %>% 
  distinct(patient_id,index_date)

# filter out 2023
index_dates <- index_dates %>% 
  filter(index_date<as.integer(ymd("2023-01-01"))) 

## Pull rx enroll locally for ccae/mdcr ----------------------------------------

tmp_collect_enroll <- function(source,year){
  db %>% 
    tbl(paste0("enrollment_detail_",source,"_",year)) %>% 
    filter(rx==1) %>% 
    select(patient_id,dtstart,dtend) %>% 
    collect() 
}

tmp_enroll_info <- collect_plan(db) %>% 
  filter(source!="medicaid") %>% 
  mutate(data = map2(source,year,tmp_collect_enroll))

load("/Volumes/AML/tmp_transfer/excess_abx/histo_med_rx_enroll.RData")

# assemble rx enrollment dates
rx_enroll <- bind_rows(tmp_enroll_info %>% 
                         unnest(data) %>% 
                         select(patient_id,dtstart,dtend),
                       select(med_rx_enroll,patient_id,dtstart,dtend)) %>% 
  arrange(patient_id,dtstart) %>% 
  mutate(gap = dtstart != lag(dtend)+1) %>% 
  mutate(gap = replace_na(gap,FALSE)) %>% 
  mutate(new_id = patient_id!=lag(patient_id)) %>% 
  mutate(new_id = replace_na(new_id,FALSE)) %>% 
  mutate(tmp =  ifelse(gap==TRUE | new_id == TRUE, 1L, 0L)) %>% 
  mutate(period = cumsum(tmp)) %>% 
  group_by(patient_id,period) %>% 
  summarise(dtstart = min(dtstart),
            dtend = max(dtend)) %>% 
  ungroup()

# Final index dates
index_dates <- index_dates %>% 
  inner_join(rx_enroll) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  mutate(rx_days_before = index_date-dtstart) %>% 
  filter(rx_days_before>=365) %>% 
  select(patient_id,index_date) 


#### Assemble Final Data -------------------------------------------------------

# distinct abx visit days (to use for analysis)
final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,upper_bound)) %>% 
  group_by(patient_id,name,days_since_index) %>% 
  summarise(daysupp = sum(daysupp)) %>% 
  ungroup() %>% 
  mutate(period = days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

final_abx_counts <- final_abx_visits %>% 
  count(name,week)

# remove incomplete weeks
week_exclude <- final_abx_visits %>% 
  distinct(week,period) %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)

final_abx_counts <- final_abx_counts %>% 
  anti_join(week_exclude)

final_abx_visits <- final_abx_visits %>% 
  anti_join(week_exclude)

tmp_out <- final_abx_visits %>% 
  distinct(patient_id,name,week) %>% 
  mutate(disease = "Histoplasmosis")

out_data <- bind_rows(out_data,tmp_out)

rm(list = ls()[!(ls() %in% c("out_data","upper_bound"))])


########################
#### Sporotrichosis ####
########################

cond_name <- "sporotrichosis"

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx") %>% 
  rename(name = `Antibiotic Name`) %>% 
  gather(key = disease, value = include, -name)

load("/Volumes/argon_home/projects/delay_diagnosis/excess_abx/data/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == "sporo") %>% 
  filter(include %in% c("Y","?")) %>% 
  filter(!(name %in% c("linezolid","cefdinir","dicloxacillin"))) %>%
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

#### Load Disease Data ---------------------------------------------------------

db <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

abx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  filter(ndcnum %in% ndc_codes) %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dates <- index_dates %>% 
  filter(time_before_index>=upper_bound) %>% 
  distinct(patient_id,index_date)

# filter out 2023
index_dates <- index_dates %>% 
  filter(index_date<as.integer(ymd("2023-01-01"))) 

## Pull rx enroll locally for ccae/mdcr ----------------------------------------

tmp_collect_enroll <- function(source,year){
  db %>% 
    tbl(paste0("enrollment_detail_",source,"_",year)) %>% 
    filter(rx==1) %>% 
    select(patient_id,dtstart,dtend) %>% 
    collect() 
}

tmp_enroll_info <- collect_plan(db) %>% 
  filter(source!="medicaid") %>% 
  mutate(data = map2(source,year,tmp_collect_enroll))

load("/Volumes/AML/tmp_transfer/excess_abx/sporo_med_rx_enroll.RData")

# assemble rx enrollment dates
rx_enroll <- bind_rows(tmp_enroll_info %>% 
                         unnest(data) %>% 
                         select(patient_id,dtstart,dtend),
                       select(med_rx_enroll,patient_id,dtstart,dtend)) %>% 
  arrange(patient_id,dtstart) %>% 
  mutate(gap = dtstart != lag(dtend)+1) %>% 
  mutate(gap = replace_na(gap,FALSE)) %>% 
  mutate(new_id = patient_id!=lag(patient_id)) %>% 
  mutate(new_id = replace_na(new_id,FALSE)) %>% 
  mutate(tmp =  ifelse(gap==TRUE | new_id == TRUE, 1L, 0L)) %>% 
  mutate(period = cumsum(tmp)) %>% 
  group_by(patient_id,period) %>% 
  summarise(dtstart = min(dtstart),
            dtend = max(dtend)) %>% 
  ungroup()

# Final index dates
index_dates <- index_dates %>% 
  inner_join(rx_enroll) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  mutate(rx_days_before = index_date-dtstart) %>% 
  filter(rx_days_before>=365) %>% 
  select(patient_id,index_date) 


#### Assemble Final Data -------------------------------------------------------

# distinct abx visit days (to use for analysis)
final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,upper_bound)) %>% 
  group_by(patient_id,name,days_since_index) %>% 
  summarise(daysupp = sum(daysupp)) %>% 
  ungroup() %>% 
  mutate(period = days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

final_abx_counts <- final_abx_visits %>% 
  count(name,week)

# remove incomplete weeks
week_exclude <- final_abx_visits %>% 
  distinct(week,period) %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)

final_abx_counts <- final_abx_counts %>% 
  anti_join(week_exclude)

final_abx_visits <- final_abx_visits %>% 
  anti_join(week_exclude)

tmp_out <- final_abx_visits %>% 
  distinct(patient_id,name,week) %>% 
  mutate(disease = "Sporotrichosis")

out_data <- bind_rows(out_data,tmp_out)


rm(list = ls()[!(ls()=="out_data")])


###############
#### Plots ####
###############


out_data %>% 
  count(week,disease) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  facet_wrap(~disease, scales = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Since Disease Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/abx_before_after.pdf")
