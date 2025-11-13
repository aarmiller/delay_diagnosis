
library(tidyverse)
library(lubridate)
library(smallDB)
library(codeBuildr)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "blasto"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 63 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

### Connect to db
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

### Load index cases
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),"index_cases.RData"))

num_cores <- 50
cp_selected <- delay_params$cp - 1 # minus 1 as the risk factors are for within delay window

#update demo1 and demo2
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
demo2 <- demo2 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")

rural_ids <- rural_visits %>% inner_join(index_cases %>% distinct(patient_id)) %>% distinct(patient_id)

# update all_dx_visits
load(paste0(delay_base_path,"all_dx_visits.RData"))
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

# update sim_obs
sim_obs <- sim_obs %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<0)

# update time map
load(paste0(delay_base_path,"delay_tm.RData"))

problem_patient_ids <- tm %>% filter(days_since_index == 0) %>% distinct(patient_id) %>%
  anti_join(tm %>% distinct(patient_id),.)

tm <- tm %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

#Load in the setting types
load(paste0(delay_base_path,"visit_info.RData"))

# update tm_stdplac
tm_stdprov <- tm_stdprov %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>%
  filter(days_since_index<=0)

load(paste0(sim_out_path,"/sim_res_ssd.RData"))
sim_res_ssd <- sim_res_ssd %>% mutate(trial = row_number()) %>% 
  select(trial, boot_trial, res) %>% 
  unnest(res) %>% 
  select(-boot_id)

load(paste0(sim_out_path,"sim_obs_reduced.RData"))
sim_obs_reduced <- sim_obs_reduced %>% 
  mutate(days_since_index = -period) %>% 
  select(-period)

sim_res_sim_obs <- sim_res_ssd %>% 
  inner_join(sim_obs_reduced, by = "obs") 

load(paste0(sim_out_path,"boot_data.RData"))
boot_data <- boot_data %>% select(boot_trial, boot_sample)

n_trials <- nrow(distinct(sim_res_ssd,trial))

# Demo data
reg_demo <- demo1 %>% 
  mutate(female=(sex==2),
         age = index_year-dobyr,
         stdrace = as.numeric(stdrace)) %>% 
  left_join(tibble(stdrace = c(0,1,2,4,9),
                   race = c("Missing/Unknown","White","Black","Hispanic","Other")),
            by = "stdrace") %>% 
  mutate(race = fct_relevel(race,"White"))

reg_demo <- reg_demo %>% 
  left_join(demo2 %>% 
              filter(index_date<=dtend & index_date>=dtstart) %>% 
              mutate(msa_new = msa %in% c("0","")) %>% 
              mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
              mutate(source = as.factor(source)) %>% 
              distinct(patient_id,source,msa=msa_new),
            by = "patient_id")

age_cats <- c(-1,17,34,44,54,64,130)

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,age,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))


rx_visits <- con %>% tbl("all_rx_visits") %>% 
  filter(patient_id %in% local(unique(index_cases$patient_id))) %>% 
  collect()

rx_visits <- rx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-(delay_params$upper_bound)))

gc()


tm <- tm %>% rename(admdate = svcdate)

## Prepare abx indicators --------------------------------------------------

# all_abx_codes <- read_csv(paste0(delay_params$out_path,"all_abx_codes.csv"))
all_abx_codes <- read_csv(paste0(stringr::str_replace(delay_params$out_path,
                                     paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
                "all_abx_codes.csv"))
all_abx_codes <- all_abx_codes %>%
  filter(!drug_name %in% tolower(c("Nitrofurantoin", "Metronidazole", "Trimethoprim",
                                   "Cefadroxil", "Rifampin", "Methenamine", "Demeclocycline",
                                   "Rifaximin", "Dapsone", "Gemifloxacin", 
                                   "Norfloxacin", "Fosfomycin", "Ofloxacin",
                                   "Carbenicillin", "Gatifloxacin", "Cefmetazole",
                                   "Cillin", "Fidaxomicin", "Sulfacetamide", "Sulfadiazine",
                                   "Sulfadimethoxine", "Sulfadoxine"))) %>% 
  distinct(ndcnum)

abx_rx_ids <- rx_visits %>% 
  filter(ndcnum %in% all_abx_codes$ndcnum) %>% 
  filter(date<index_date & date>=(index_date-cp_selected)) %>%
  distinct(patient_id) %>% 
  mutate(resp_antibiotic_drugs_window=1L)

reg_demo <- reg_demo %>%
  left_join(abx_rx_ids) %>% 
  mutate_at(vars(resp_antibiotic_drugs_window),~replace_na(.,0L))

## Prepare inhaler indicators --------------------------------------------------

tmp <- read_csv("/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/data/all_inhaler_codes.csv")

inhale_rx_ids <- rx_visits %>% 
  filter(ndcnum %in% tmp$NDCNUM) %>% 
  filter(date<index_date & date>=(index_date-cp_selected)) %>%
  distinct(patient_id) %>% 
  mutate(inhalers_window=1L)

reg_demo <- reg_demo %>%
  left_join(inhale_rx_ids) %>% 
  mutate_at(vars(inhalers_window),~replace_na(.,0L))

## Prepare prior asthma and COPD indicators ------------------------------------

#load asthma/copd codes
old_ssds <- read_rds(paste0("/Shared/Statepi_Diagnosis/params/ssd_codes/condition_specific/ssd_tb.RDS"))
old_ssds <- old_ssds[c("asthma", "copd")]
asthma_COPD_codes <- tibble(dx = old_ssds$asthma$icd9_codes, dx_ver= 9L, cond = "Asthma") %>% 
  bind_rows(tibble(dx = old_ssds$asthma$icd10_codes, dx_ver= 10L, cond = "Asthma")) %>% 
  bind_rows(tibble(dx = old_ssds$copd$icd9_codes, dx_ver= 9L, cond = "COPD")) %>% 
  bind_rows(tibble(dx = old_ssds$copd$icd10_codes, dx_ver= 10L, cond = "COPD")) %>% 
  distinct()

dx_visits <- con %>% tbl("all_dx_visits") %>% collect()

dx_visits <- dx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)

#Asthma
asthma <- dx_visits %>% 
  inner_join(asthma_COPD_codes %>% filter( cond == "Asthma") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(asthma_prior_cp = 1L) 

#COPD
copd <- dx_visits %>% 
  inner_join(asthma_COPD_codes %>% filter( cond == "COPD") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(copd_prior_cp = 1L) 

reg_demo <- reg_demo %>%
  left_join(asthma) %>% 
  left_join(copd) %>% 
  mutate_at(vars(asthma_prior_cp, copd_prior_cp),~replace_na(.,0L))

## Prepare prior chest imaging -------------------------------------------------
chest_ct <-c("71260", "71250", "71270")
chest_xray <- c("71010", "71015", "71020", "71021", "71022", "71023",
                "71030", "71034", "71035", "71101", "71111", "71120",
                "71045", "71046", "71047", "71048")

proc_visits <- con %>% tbl("all_proc_visits") %>% collect()

proc_visits <- proc_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)

chest_ct_visits <- proc_visits %>% 
  filter(proc %in% chest_ct)%>% 
  distinct(patient_id) %>% mutate(chest_ct_prior_cp = 1)

chest_xray_visits <- proc_visits %>% 
  filter(proc %in% chest_xray) %>% 
  distinct(patient_id) %>% mutate(chest_xray_prior_cp = 1)

reg_demo <- reg_demo %>%
  left_join(chest_ct_visits) %>% 
  left_join(chest_xray_visits) %>% 
  mutate_at(vars(chest_ct_prior_cp, chest_xray_prior_cp),~replace_na(.,0L))

## Prepare immunosuppression, hiv, and diabetes indicators ------------------------------------

#code
# IMMUNOSUPRESSION (lump all codes together).
# 
# ICD-10-CM   Z79.5 Long term (current) use of steroids
# ICD-9-CM    V58.65 Long-term (current) use of steroids
# 
# ICD-9-CM   V87.46 Personal history of immunosuppressive therapy
# ICD-10-CM   Z92.25 Personal history of immunosuppression therapy
# 
# ICD-9-CM    no code for Long term (current) use of immunomodulators and immunosuppressants
# ICD-10-CM   Z79.6 Long term (current) use of immunomodulators and immunosuppressants
# 
# ICD-9-CM    042 Human immunodeficiency virus [HIV] disease
# ICD-10-CM   B20 Human immunodeficiency virus [HIV] disease
# 
# 
# DIABETES:
#   
#   ICD-9-CM    250 Diabetes mellitus
# 
# ICD-10-CM   E10 Type 1 diabetes mellitus
# ICD-10-CM   E11 Type 2 diabetes mellitus

immuno_dia_hiv_codes <- tibble(dx_supp = c("Z795", "V5865",
                                      "V8746", "Z9225",
                                      "Z796", "Z7960", "Z7961", "Z7962", "Z79620", "Z79621", "Z79622", "Z79623", "Z79624",
                                      "042", "B20",
                                      "250", "E10", "E11" ),
                               dx_ver = c(10,9,9,10,rep(10, 9),9,10,9,10,10),
                               cond = c(rep("Immunosuppression", 13),
                                        rep("HIV", 2),
                                        rep("Diabetes", 3))) %>% 
  mutate(dx = map2(dx_supp, dx_ver, ~{
    if(.y == 9){
      children_safe(.x)
    } else {
      temp <- children10(.x)
      if(length(temp)==0){
        .x
      } else {
        temp
      }
    }})) %>% 
  unnest(dx) %>% 
  select(-dx_supp) %>% 
  distinct()

dx_visits <- con %>% tbl("all_dx_visits") %>% collect()

dx_visits <- dx_visits %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(days_since_index< (-cp_selected) & days_since_index>=-delay_params$upper_bound)

# Immunosuppression (40)
immunosuppression <- dx_visits %>% 
  inner_join(immuno_dia_hiv_codes %>% filter( cond == "Immunosuppression") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(immunosuppression_prior_cp = 1L) 

# HIV (40) only one in common with Immunosuppression 
hiv <- dx_visits %>% 
  inner_join(immuno_dia_hiv_codes %>% filter( cond == "HIV") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(hiv_prior_cp = 1L) 

# Diabetes (670) 
diabetes <- dx_visits %>% 
  inner_join(immuno_dia_hiv_codes %>% filter( cond == "Diabetes") %>% select(-cond)) %>% 
  distinct(patient_id) %>% mutate(diabetes_prior_cp = 1L) 

reg_demo <- reg_demo %>%
  left_join(immunosuppression) %>% 
  left_join(hiv) %>% 
  left_join(diabetes) %>% 
  mutate_at(vars(immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp),~replace_na(.,0L))

## Prepare IPD 

#IPD
ipd_codes <- codeBuildr::load_disease_codes("ipd")
ipd_codes <- tibble(dx_ver = 9, dx = ipd_codes$ipd$icd9_codes) %>% 
  bind_rows(tibble(dx_ver = 10, dx = ipd_codes$ipd$icd10_codes))

IPD <- dx_visits %>% inner_join(ipd_codes) %>%   
  distinct(patient_id) %>% mutate(IPD_prior_cp = 1L) 

reg_demo <- reg_demo %>%
  left_join(IPD) %>% 
  mutate_at(vars(IPD_prior_cp),~replace_na(.,0L))

## Prepare prior geography  ----------------------------------------------------
# source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
# load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024

source(paste0(stringr::str_replace(delay_params$out_path,
                                   paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
              "get_enroll_detail_fun.R"))
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
            "egeoloc_labels.RData"))

enroll_collapsed_temp <- gather_collapse_enrollment(enrolid_list = reg_demo %>% distinct(patient_id) %>% .$patient_id,
                                                    vars = "egeoloc",
                                                    db_path =  paste0(delay_params$small_db_path,"blasto.db"),
                                                    num_cores=10,
                                                    collect_tab = collect_table(year = 1:22))

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322071/
enroll_collapsed_temp2 <- enroll_collapsed_temp %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  inner_join(select(demo1,patient_id,index_date)) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  distinct() 
# only 3,567 out of 3,825 have location information

# top2_high_inc_state_baddley baddely states with incidence >= 1.38

location_ind <- enroll_collapsed_temp2 %>%  
  mutate(top2_high_inc_state_baddley = ifelse(state_abb %in% c("ND", "MN",
                                                               "WI", "IL",
                                                               "TN", "AL",
                                                               "MS", "LA"), 1L, ifelse(is.na(state_abb), NA, 0))) %>% 
  select(patient_id, location:state_abb, top2_high_inc_state_baddley) %>% 
  distinct()
# location_ind %>% filter(top2_high_inc_state_baddley == 1) %>% distinct(state_abb) #check coding
# 3,387 of the  3,567 with location info have non missing state

reg_demo <- reg_demo %>%
  left_join(location_ind) 

  # need to keep NAs, either medicaid or those with no state info
  # mutate(top2_high_inc_state_baddley = ifelse(source!="medicaid" & is.na(top2_high_inc_state_baddley), 0L, top2_high_inc_state_baddley))
  # mutate_at(vars(top2_high_inc_state_baddley),~replace_na(.,0L)) # need to keep NAs for medicaid


## Top 5 truven inc states
load(paste0(out_path,"truven_inc_by_egeoloc_year.RData")) 

top5_inc_states <- ave_truven_blasto_inc_state %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  filter(!is.na(state_abb)) %>% 
  arrange(desc(ave_inc_per_100000)) %>% 
  filter(row_number()<=5) %>% 
  .$state_abb

location_ind2 <- enroll_collapsed_temp2 %>%  
  mutate(top5_truven_inc_states = ifelse(state_abb %in% top5_inc_states, 1L, ifelse(is.na(state_abb), NA, 0))) %>% 
  select(patient_id, location:state_abb, top5_truven_inc_states) %>% 
  distinct()

reg_demo <- reg_demo %>%
  left_join(location_ind2) 


## SAVE COUNTS
index_state_counts <- reg_demo %>% 
  group_by(location, state_name, state_abb) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  mutate(percent = n/sum(n)*100) %>% 
  arrange(desc(n)) 
write_csv(index_state_counts, file = paste0(out_path,"index_counts_by_state.csv"))  

index_top2_high_inc_state_baddley_count <- reg_demo %>% 
  count(top2_high_inc_state_baddley) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  arrange(desc(n)) 
write_csv(index_top2_high_inc_state_baddley_count, file = paste0(out_path,"index_top2_high_inc_state_baddley_count.csv"))  

index_top5_truven_inc_states <- reg_demo %>% 
  count(top5_truven_inc_states) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  arrange(desc(n)) 
write_csv(index_top5_truven_inc_states, file = paste0(out_path,"index_top5_truven_inc_states.csv"))  

index_WI_MN_count <- reg_demo %>% 
  mutate(WI_MN = ifelse(is.na(state_abb), NA, (state_abb %in% c("WI", "MN")) *1.0)) %>% 
  count(WI_MN) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  arrange(desc(n)) 
write_csv(index_WI_MN_count, file = paste0(out_path,"index_WI_MN_count.csv"))  

west_states <- c("WA", "OR", "CA",
                "ID", "NV", "AZ",
                "MT", "WY", "CO",
                "UT", "NM")

index_west_states_count <- reg_demo %>% 
  mutate(west_states = ifelse(is.na(state_abb), NA, (state_abb %in% west_states)*1.0)) %>%
  count(west_states) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  arrange(desc(n)) 
write_csv(index_west_states_count, file = paste0(out_path,"index_west_states_count.csv"))  

## Add in ID consults ----------------------------------------------------------
ID_consult_vis <- tm_stdprov %>% 
  filter(patient_id %in% local(index_cases$patient_id)) %>%
  filter(stdprov %in%  c("285", "448")) %>%
  distinct() %>%
  inner_join(index_cases %>%
               mutate(index_date = index_date + shift) %>%
               select(patient_id, index_date), by = "patient_id") %>%
  rename(admdate = svcdate) %>%
  filter(admdate<=index_date & (admdate)>=(index_date-(delay_params$upper_bound))) %>%
  distinct(patient_id, admdate)

tm <- tm  %>%
  left_join(ID_consult_vis %>% mutate(ID_consult = 1L), by = c("patient_id", "admdate")) %>% 
  mutate_at(vars(ID_consult),~replace_na(.,0L))


############################
#### Prepare Visit Info ####
############################

# total number of combinations: \sum_i=1^k choose(n, k) + 1 for no cat
# so here it is sum(choose(4, 1:4)) +1

setting_labels <- expand.grid(outpatient = c("Outpatient", NA),
                    inpatient = c("Inpatient", NA),
                    ed = c("ED", NA),
                    obs_stay = c("Observational Stay", NA)) %>% 
  unite(., col = setting_label, sep = " + ", na.rm = T, remove =F) %>% 
  mutate(across(outpatient:obs_stay, ~ifelse(is.na(.), F, T))) %>% 
  mutate(setting_label = ifelse(setting_label == "", "Not any",
                                ifelse(setting_label == "Inpatient", "Inpatient Only",
                                       ifelse(setting_label == "Outpatient", "Outpatient Only",
                                              ifelse(setting_label == "ED", "ED Only",
                                                     ifelse(setting_label == "Observational Stay", "Observational Stay Only",
                                                            setting_label))))))
setting_labels <- setting_labels[order(rowSums(setting_labels[,2:5])), ]
setting_labels <- setting_labels %>% mutate(setting_label = factor(setting_label, 
                                 levels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")]),
                                 labels = c(setting_labels$setting_label[which(setting_labels$setting_label=="Outpatient Only")],
                                            setting_labels$setting_label[-which(setting_labels$setting_label %in% c("Not any", "Outpatient Only"))],
                                            setting_labels$setting_label[which(setting_labels$setting_label=="Not any")])))

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(patient_id,outpatient,ed,obs_stay,inpatient,admdate, ID_consult)) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, ID_consult) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(obs,patient_id,outpatient,ed,obs_stay,inpatient,admdate, ID_consult)) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(obs,patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month, ID_consult) %>% 
  inner_join(setting_labels) 

### Prep weekend visit info ------------

# add weekend and demo to sim data
sim_res_ssd <- sim_res_ssd %>% 
  inner_join(obs_locations)


full_reg_data <- nest(sim_res_ssd, .by = c("trial", "boot_trial"))
full_reg_data <- full_reg_data %>% 
  inner_join(boot_data %>% 
               mutate(boot_sample = map(boot_sample, ~inner_join(., index_locations))))

# all.equal(full_reg_data %>% slice(100) %>% select(boot_sample) %>% unnest(boot_sample),
#           full_reg_data %>% slice(1000) %>% select(boot_sample) %>% unnest(boot_sample))

print("dataset ready")

###########################
#### Regression Models ####
###########################

# exclude year 2001 due to sparsity


#### Missed opportunities All (NO SETTING RISK FACTOR) --------------------------------------------------

get_miss_res <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, weekend, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  fit <- glm(miss ~ .,
             family = "binomial", data=reg_data)
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit))
  
}

miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_miss_res(x)}, 
                                   mc.cores = num_cores)


miss_opp_res <- bind_rows(miss_opp_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

# rm(cluster)
gc()

print("mod 1a finished")

#### Missed opportunities All ID_consult and top2_high_inc_state_baddley INTERACTION (NO SETTING RISK FACTOR) --------------------------------------------------
get_miss_res_interaction <- function(trial_val){
  
  # trial_val <- 1
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year>2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, weekend, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  rhs_flm <- paste0(paste0(colnames(reg_data %>% select(-miss)), collapse = " + "), " + ", "top2_high_inc_state_baddley:ID_consult")
  flm <- as.formula(paste0("miss ~ ", rhs_flm))
  
  fit <- glm(flm,
             family = "binomial", data=reg_data)
  
  coefs <- coef(fit)
  effect_state_by_ID <- tibble(term = c("Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - ID_consult",
                                        "Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - No ID_consult",
                                        "Est ID_consult vs no ID_consult - top2_high_inc_state_baddley",
                                        "Est ID_consult vs no ID_consult - Non-top2_high_inc_state_baddley"
  ),
  estimate = c(coefs["top2_high_inc_state_baddley"] + coefs["top2_high_inc_state_baddley:ID_consult"],
               coefs["top2_high_inc_state_baddley"],
               coefs["ID_consult"] + coefs["top2_high_inc_state_baddley:ID_consult"],
               coefs["ID_consult"]))
  
  tibble(term = names(coef(fit)),
         estimate = coef(fit)) %>% 
    bind_rows(effect_state_by_ID)
  
}

miss_opp_res_interaction <- parallel::mclapply(1:max(full_reg_data$trial),
                                               function(x){get_miss_res_interaction(x)}, 
                                               mc.cores = num_cores)

miss_opp_res_interaction <- bind_rows(miss_opp_res_interaction) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

# rm(cluster)
gc()

print("mod 1b finished")


# #### Missed opportunities All --------------------------------------------------
# library(logistf)
# get_miss_res <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(year>2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(miss, setting_label, age_cat, female, rural, source, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
#     drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
#   
#   fit <- logistf(miss ~ .,
#                  family = "binomial", data=reg_data,
#                  control = logistf.control( maxit = 1000),
#                  plcontrol = logistpl.control( maxit = 1000))
#   
#   tibble(term = names(coef(fit)),
#          estimate = coef(fit))
#   
# }
# 
# miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
#                                    function(x){get_miss_res(x)}, 
#                                    mc.cores = num_cores)
# 
# 
# miss_opp_res <- bind_rows(miss_opp_res) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
#             low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
#             high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))
# 
# # rm(cluster)
# gc()
# 
# print("mod 1a finished")
# 
# #### Missed opportunities All ID_consult and top2_high_inc_state_baddley INTERACTION --------------------------------------------------
# get_miss_res_interaction <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(year>2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(miss, setting_label, age_cat, female, rural, source, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
#     drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
#   
#   rhs_flm <- paste0(paste0(colnames(reg_data %>% select(-miss)), collapse = " + "), " + ", "top2_high_inc_state_baddley:ID_consult")
#   flm <- as.formula(paste0("miss ~ ", rhs_flm))
# 
#   fit <- logistf(flm,
#                  family = "binomial", data=reg_data,
#                  control = logistf.control( maxit = 1000),
#                  plcontrol = logistpl.control( maxit = 1000))
#   
#   coefs <- coef(fit)
#   effect_state_by_ID <- tibble(term = c("Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - ID_consult",
#                                         "Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - No ID_consult",
#                                         "Est ID_consult vs no ID_consult - top2_high_inc_state_baddley",
#                                         "Est ID_consult vs no ID_consult - Non-top2_high_inc_state_baddley"
#                                         ),
#                                estimate = c(coefs["top2_high_inc_state_baddley"] + coefs["top2_high_inc_state_baddley:ID_consult"],
#                                             coefs["top2_high_inc_state_baddley"],
#                                             coefs["ID_consult"] + coefs["top2_high_inc_state_baddley:ID_consult"],
#                                             coefs["ID_consult"]))
#   
#   tibble(term = names(coef(fit)),
#          estimate = coef(fit)) %>% 
#     bind_rows(effect_state_by_ID)
#   
# }
# 
# miss_opp_res_interaction <- parallel::mclapply(1:max(full_reg_data$trial),
#                                                function(x){get_miss_res_interaction(x)}, 
#                                                mc.cores = num_cores)
# 
# miss_opp_res_interaction <- bind_rows(miss_opp_res_interaction) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
#             low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
#             high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))
# 
# # rm(cluster)
# gc()
# 
# print("mod 1b finished")
# 
# #### Missed opportunities All -inpatient only ind --------------------------------------------------
# 
# get_miss_res_inpatient_ind <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(year > 2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(miss, inpatient, age_cat, female, rural, source, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
#     drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
#   
#   
#   fit <- glm(miss ~ .,
#              family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# miss_opp_res_inpatient_ind <- parallel::mclapply(1:max(full_reg_data$trial),
#                                                   function(x){get_miss_res_inpatient_ind(x)}, 
#                                                   mc.cores = num_cores)
# 
# 
# miss_opp_res_inpatient_ind <- bind_rows(miss_opp_res_inpatient_ind) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
#             low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
#             high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))
# 
# # rm(cluster)
# gc()
# 
# print("mod 2a finished")
# 
# 
# #### Missed opportunities All -inpatient only ind ID_consult and top2_high_inc_state_baddley INTERACTION --------------------------------------------------
# 
# get_miss_res_inpatient_ind_interaction <- function(trial_val){
#   
#   # trial_val <- 1
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE)) %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(year > 2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(miss, inpatient, age_cat, female, rural, source, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, top2_high_inc_state_baddley, ID_consult) %>% 
#     drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
#   
#   rhs_flm <- paste0(paste0(colnames(reg_data %>% select(-miss)), collapse = " + "), " + ", "top2_high_inc_state_baddley:ID_consult")
#   flm <- as.formula(paste0("miss ~ ", rhs_flm))
#   
#   fit <- glm(flm,
#              family = "binomial", data=reg_data)
#   
#   coefs <- coef(fit)
#   effect_state_by_ID <- tibble(term = c("Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - ID_consult",
#                                         "Est top2_high_inc_state_baddley vs Non-top2_high_inc_state_baddley - No ID_consult",
#                                         "Est ID_consult vs no ID_consult - top2_high_inc_state_baddley",
#                                         "Est ID_consult vs no ID_consult - Non-top2_high_inc_state_baddley"
#   ),
#   estimate = c(coefs["top2_high_inc_state_baddley"] + coefs["top2_high_inc_state_baddley:ID_consult"],
#                coefs["top2_high_inc_state_baddley"],
#                coefs["ID_consult"] + coefs["top2_high_inc_state_baddley:ID_consult"],
#                coefs["ID_consult"]))
#   
#   tibble(term = names(coef(fit)),
#          estimate = coef(fit)) %>% 
#     bind_rows(effect_state_by_ID)
#   
# }
# 
# miss_opp_res_inpatient_ind_interaction <- parallel::mclapply(1:max(full_reg_data$trial),
#                                                  function(x){get_miss_res_inpatient_ind_interaction(x)}, 
#                                                  mc.cores = num_cores)
# 
# 
# miss_opp_res_inpatient_ind_interaction <- bind_rows(miss_opp_res_inpatient_ind_interaction) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
#             low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
#             high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))
# 
# # rm(cluster)
# gc()
# 
# print("mod 2b finished")

# #### Missed opportunities Medicaid ---------------------------------------------
# 
# get_miss_res_med <- function(trial_val){
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE))  %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#     filter(year > 2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     filter(source == "medicaid") %>% 
#     select(miss, setting_label, age_cat, female, rural, race, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp)
#   
#   
#   fit <- glm(miss~.,
#              family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# 
# miss_opp_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
#                                             function(x){get_miss_res_med(x)}, 
#                                             mc.cores = num_cores)
# 
# miss_opp_res_med <- bind_rows(miss_opp_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# gc()
# 
# print("mod 3 finished")
# 
# #### Missed opportunities Medicaid - Inpatient_only ---------------------------------------------
# 
# get_miss_res_med_inpatient_ind <- function(trial_val){
#   
#   tmp1 <- full_reg_data %>% 
#     filter(trial==trial_val)
#   
#   reg_data <- bind_rows(tmp1 %>% select(data) %>% 
#                           unnest(data) %>% 
#                           mutate(miss=TRUE),
#                         tmp1 %>% select(boot_sample) %>% 
#                           unnest(boot_sample) %>% 
#                           mutate(miss=FALSE))  %>% 
#     inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
#      filter(year > 2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     filter(source == "medicaid") %>% 
#     select(miss, inpatient, age_cat, female, rural, race, weekend, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp)
#   
#   
#   fit <- glm(miss ~ .,
#              family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# 
# miss_opp_res_med_inpatient_ind  <- parallel::mclapply(1:max(full_reg_data$trial),
#                                        function(x){get_miss_res_med_inpatient_ind(x)}, 
#                                        mc.cores = num_cores)
# 
# miss_opp_res_med_inpatient_ind <- bind_rows(miss_opp_res_med_inpatient_ind) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# gc()
# 
# print("mod 4 finished")

#### Delay duration All --------------------------------------------------------

# compute duration by simulation (note this will be used in the next step)
full_reg_data_dur <- full_reg_data %>% 
  mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")) %>% 
                      group_by(patient_id) %>% 
                      summarise(duration = -min(days_since_index)) %>% 
                      ungroup()))

get_dur_res <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(duration, age_cat, female, rural, source, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, #top2_high_inc_state_baddley,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
  fit <- glm(duration~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                       function(x){get_dur_res(x)}, 
                                       mc.cores = num_cores)


miss_dur_res <- bind_rows(miss_dur_res) %>% 
  group_by(term) %>% 
  summarise(est = median(estimate, na.rm = T),
            low = quantile(estimate, probs = c(0.025), na.rm = T),
            high = quantile(estimate, probs = c(0.975), na.rm = T),
            low_90 = quantile(estimate, probs = c(0.05), na.rm = T),
            high_90 = quantile(estimate, probs = c(0.95), na.rm = T) )


gc()

print("mod 3 finished")


# log normal model
get_dur_res_log_normal <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    mutate(log_duration = log(duration)) %>% 
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(log_duration, age_cat, female, rural, source, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, #top2_high_inc_state_baddley,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
  fit <- glm(log_duration ~., family = "gaussian", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_log_normal <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_dur_res_log_normal(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_log_normal <- bind_rows(miss_dur_res_log_normal) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))


gc()

print("mod 4 finished")

# weibul model
library(survival)
get_dur_res_weibull <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    mutate(delta = 1) %>% 
    select(duration, delta, age_cat, female, rural, source, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, #top2_high_inc_state_baddley,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() %>%  # have to do at the end as obs not in bootsample and boot_id not in data
    mutate(source = as.factor(as.character(source))) # for some reason survreg will not drop the missing factor level medicaid
                                       # and when your run summary fit the NA coef will be given a value of 0
  
  fit <- survreg(Surv(duration, delta) ~., dist = "weibull", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_dur_res_weibull <- parallel::mclapply(1:max(full_reg_data$trial),
                                   function(x){get_dur_res_weibull(x)}, 
                                   mc.cores = num_cores)


miss_dur_res_weibull <- bind_rows(miss_dur_res_weibull) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

gc()

print("mod 5 finished")

# #### Delay duration Medicaid ---------------------------------------------------
# 
# get_dur_res_med <- function(trial_val){
#   
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     select(patient_id) %>%     
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L)) %>% 
#     filter(source == "medicaid") %>% 
#     filter(year > 2001) %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(duration, age_cat, female, rural, race, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp,
#            resp_antibiotic_drugs_window, inhalers_window)
#   
#   
#   fit <- glm(duration~., family = "gaussian", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# 
# miss_dur_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
#                                    function(x){get_dur_res_med(x)}, 
#                                    mc.cores = num_cores)
# 
# 
# miss_dur_res_med <- bind_rows(miss_dur_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(estimate, na.rm = T),
#             low = quantile(estimate ,probs = c(0.025), na.rm = T),
#             high = quantile(estimate ,probs = c(0.975), na.rm = T))
# 
# gc()
# 
# print("mod 6 finished")


#### Delay Patient All ---------------------------------------------------------


get_delay_pat_res <- function(trial_val){
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>%     
    inner_join(reg_demo, by = "patient_id") %>% 
    left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
    mutate(duration = replace_na(duration,0L)) %>%  
    mutate(miss = duration>0) %>% 
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, age_cat, female, rural, source, year, month, immunosuppression_prior_cp, hiv_prior_cp, diabetes_prior_cp,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, #top2_high_inc_state_baddley,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  
  
  fit <- glm(miss~., family = "binomial", data=reg_data)
  
  broom::tidy(fit)
  
  
}


miss_delay_pat_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                       function(x){get_delay_pat_res(x)}, 
                                       mc.cores = num_cores)

miss_delay_pat_res <- bind_rows(miss_delay_pat_res) %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate),probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T))

gc()

print("mod 6 finished")

# #### Delay Patient Medicaid ----------------------------------------------------
# 
# 
# get_delay_pat_res_med <- function(trial_val){
# 
#   tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
#   
#   reg_data <- tmp1 %>% 
#     select(boot_sample) %>% 
#     unnest(boot_sample) %>% 
#     select(patient_id) %>%     
#     inner_join(reg_demo, by = "patient_id") %>% 
#     left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>% 
#     mutate(duration = replace_na(duration,0L)) %>%  
#     mutate(miss = duration>0) %>% 
#     filter(year > 2001) %>% 
#     filter(source == "medicaid") %>% 
#     mutate(year = as.factor(year),
#            month = factor(month, levels = 1:12)) %>% 
#     select(miss, age_cat, female, rural, race, year, month,
#            asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp,
#            resp_antibiotic_drugs_window, inhalers_window)
#   
#   
#   fit <- glm(miss~., family = "binomial", data=reg_data)
#   
#   broom::tidy(fit)
#   
#   
# }
# 
# miss_delay_pat_res_med <- parallel::mclapply(1:max(full_reg_data$trial),
#                                          function(x){get_delay_pat_res_med(x)}, 
#                                          mc.cores = num_cores)
# 
# miss_delay_pat_res_med <- bind_rows(miss_delay_pat_res_med) %>% 
#   group_by(term) %>% 
#   summarise(est = median(exp(estimate), na.rm = T),
#             low = quantile(exp(estimate),probs = c(0.025), na.rm = T),
#             high = quantile(exp(estimate),probs = c(0.975), na.rm = T))
# 
# 
# gc()
# 
# print("mod 8 finished")

########################
##### Save Results #####
########################

ssd_miss_risk_models <- list(miss_opp_res = miss_opp_res,
                             miss_opp_res_interaction = miss_opp_res_interaction,
                             # miss_opp_res_inpatient_ind = miss_opp_res_inpatient_ind,
                             # miss_opp_res_inpatient_ind_interaction = miss_opp_res_inpatient_ind_interaction,
                             # miss_opp_res_med = miss_opp_res_med,
                             # miss_opp_res_med_inpatient_ind = miss_opp_res_med_inpatient_ind,
                             miss_dur_res = miss_dur_res,
                             miss_dur_res_log_normal = miss_dur_res_log_normal,
                             miss_dur_res_weibull = miss_dur_res_weibull,
                             # miss_dur_res_med = miss_dur_res_med,
                             miss_delay_pat_res = miss_delay_pat_res
                             # miss_delay_pat_res_med = miss_delay_pat_res_med
                             )


save(ssd_miss_risk_models, index_locations, file = paste0(out_path,"ssd_miss_risk_models.RData"))
# load(paste0(out_path,"ssd_miss_risk_models.RData"))


rmarkdown::render(input = paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/", proj_name, "/risk_model_report.Rmd"),
                  params = list(cond = tools::toTitleCase(proj_name)),
                  output_dir = out_path,
                  output_file = paste0(proj_name, "_risk_model_report_", lubridate::today() %>% format('%m-%d-%Y')))

