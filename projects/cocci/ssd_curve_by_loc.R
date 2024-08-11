rm(list=ls())
library(tidyverse)
library(bit64)

# load final delay param
cond_name <- "cocci"
load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[cond_name]]
out_path <-paste0(delay_params$out_path,"risk_models/") 
if (!dir.exists(delay_params$out_path)) {
  dir.create(delay_params$out_path)
}

## Build the primary index_dx_visits -------------------------------------------
# Connect to cocci DB
db <- src_sqlite(paste0(delay_params$small_db_path, cond_name, ".db"))
load(paste0(delay_params$out_path,"index_cases.RData"))

## Prepare prior geography  ----------------------------------------------------
delay_base_path <- paste0(delay_params$base_path,"delay_results/")
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
demo2 <- demo2 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")

source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024

enroll_collapsed_temp <- gather_collapse_enrollment(enrolid_list = index_cases %>% distinct(patient_id) %>% .$patient_id,
                                                    vars = "egeoloc",
                                                    db_path =  paste0(delay_params$small_db_path,"cocci.db"),
                                                    num_cores=10,
                                                    collect_tab = collect_table(year = 1:22))

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322071/
enroll_collapsed_temp2 <- enroll_collapsed_temp %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  distinct() 
# only 26,378 out of 26,905  have location information

# top2_high_inc_state_baddley baddely states with incidence >= 3.97

location_ind <- enroll_collapsed_temp2 %>%  
  mutate(top2_high_inc_state_baddley = ifelse(state_abb %in% c("ID", "WY",
                                                               "CA", "NV",
                                                               "UT", "AZ",
                                                               "NM", "ND",
                                                               "MN", "IA",
                                                               "AR", "NH",
                                                               "AK"), 1L, ifelse(is.na(state_abb), NA, 0))) %>% 
  select(patient_id, location:state_abb, top2_high_inc_state_baddley) %>% 
  distinct()
# location_ind %>% filter(top2_high_inc_state_baddley == 1) %>% distinct(state_abb) #check coding
# 26,048 of the 26,378 with location info have non missing state

loc_index_cases <- index_cases %>% select(patient_id) %>% 
  left_join(location_ind) %>% 
  count(top2_high_inc_state_baddley)


# update all_dx_visits
load(paste0(delay_base_path,"all_dx_visits.RData"))
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

SSD_codes <- codeBuildr::load_ssd_codes("cocci")
SSD_codes <- SSD_codes %>% mutate(dx_ver = ifelse(type == "icd9", 9, 10)) %>% 
  rename(dx = code) %>% 
  select(-type)

ssd_curve_w_loc <- all_dx_visits %>% 
  inner_join(SSD_codes) %>% 
  distinct(patient_id, days_since_index) %>% 
  left_join(loc_index_cases) 
  
save(ssd_curve_w_loc, file = paste0(out_path,"ssd_curve_w_loc.RData"))  
# 
# out_path <- "/Volumes/Statepi_Diagnosis/projects/cocci/risk_models/"
# load(paste0(out_path, "/ssd_curve_w_loc.RData"))
# 
# ssd_curve_w_loc %>% 
#   filter(days_since_index<0) %>% 
#   group_by(top2_high_inc_state_baddley) %>% 
#   count(days_since_index) %>% 
#   ungroup() %>% 
#   mutate(group = ifelse(is.na(top2_high_inc_state_baddley), "Missing", top2_high_inc_state_baddley)) %>% 
#   mutate(group = factor(group, levels = c("0","1","Missing"), 
#                                               labels = c("No", "Yes", "Missing"))) %>% 
#   # count(group, top2_high_inc_state_baddley)
#   bind_rows(ssd_curve_w_loc %>% 
#               filter(days_since_index<0) %>% 
#               count(days_since_index) %>% 
#               mutate(group = "Overall")) %>% 
#   ggplot(aes(y = n, x = days_since_index, col = group))+
#   geom_line()
# 
# 
# ################################################################################
# proc_visit <- tbl(db, "all_proc_visits") %>% 
#   filter(days_since_index <0 & days_since_index > -7) %>% 
#   collect(n=Inf)
# 
# temp <- proc_visit %>% 
#   distinct() %>% 
#   count(proc) %>% 
#   arrange(desc(n))
# 
# cpt_labels <- read_csv(paste0(out_path, "/common_cpt_within_7days.csv")) %>% 
#   select(proc = description, desc = status)
# 
# temp %>% 
#   left_join(cpt_labels) %>% 
#   print(n = 100)
# 
# proc_visit %>% filter(proc == "88305") %>% 
#   distinct(patient_id) %>% mutate(ind = 1) %>% 
#   left_join(loc_index_cases,.) %>% 
#   mutate(ind = ifelse(is.na(ind), 0, ind)) %>% 
#   count(top2_high_inc_state_baddley, ind) %>% 
#   group_by(top2_high_inc_state_baddley) %>% 
#   mutate(percent= n/(sum(n))*100) %>% 
#   ungroup()
#   
