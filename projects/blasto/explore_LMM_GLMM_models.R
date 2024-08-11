

library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "blasto"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

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
load(paste0(delay_params$out_path,"index_cases.RData"))

num_cores <- 50
cp_selected <- delay_params$cp

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

all_abx_codes <- read_csv(paste0(delay_params$out_path,"all_abx_codes.csv"))

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
source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024

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
  distinct_at(vars(patient_id,outpatient,ed,obs_stay,inpatient,admdate)) %>% 
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month) %>% 
  inner_join(setting_labels) 

obs_locations <- tm %>% 
  inner_join(sim_res_sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(outpatient==1 | ed==1 | obs_stay==1 | inpatient==1) %>%
  distinct_at(vars(obs,patient_id,outpatient,ed,obs_stay,inpatient,admdate)) %>%
  mutate(dow = weekdays(as_date(admdate))) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  mutate(year = year(as_date(admdate)),
         month = month(as_date(admdate))) %>% 
  select(obs,patient_id,outpatient,ed,obs_stay,inpatient,weekend,dow, year, month) %>% 
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

# compute duration by simulation (note this will be used in the next step)
full_reg_data_dur <- full_reg_data %>% 
  mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")) %>% 
                      group_by(patient_id) %>% 
                      summarise(duration = -min(days_since_index)) %>% 
                      ungroup()))

print("dataset ready")

#### Delay Patient All ---------------------------------------------------------
library(lme4)

get_delay_res_LMM <- function(trial_val){
  
  
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
    select(duration, age_cat, female, rural, source, year, month,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp,
           #top2_high_inc_state_baddley,
           state_abb,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  fit <- lme4::lmer(duration ~ age_cat + female + rural + source + year + month + 
               asthma_prior_cp + copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp + 
               resp_antibiotic_drugs_window + inhalers_window +(1|state_abb),
             data = reg_data)
  
  RE_temp <- lme4::ranef(fit)$state_abb
  FE_temp <- lme4::fixef(fit)
  
  RE <- tibble(intercept = RE_temp[,1],
               state_abb = rownames(RE_temp))
  FE <- tibble(term = names(FE_temp),
               estimate = FE_temp)
  
  return(list(RE = RE,
              FE = FE))
         
}


miss_delay_pat_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                         function(x){get_delay_res_LMM(x)}, 
                                         mc.cores = num_cores)

RE <- do.call("rbind", lapply(miss_delay_pat_res, function(x){x$RE} ))
FE <- do.call("rbind", lapply(miss_delay_pat_res, function(x){x$FE} ))


miss_delay_dur_LMM_FE <- FE %>% 
  group_by(term) %>% 
    summarise(est = median(estimate, na.rm = T),
              low = quantile(estimate, probs = c(0.025), na.rm = T),
              high = quantile(estimate, probs = c(0.975), na.rm = T),
              low_90 = quantile(estimate, probs = c(0.05), na.rm = T),
              high_90 = quantile(estimate, probs = c(0.95), na.rm = T) )

miss_delay_dur_LMM_RE <- RE %>% 
  group_by(state_abb) %>% 
  summarise(intercept_est = median(intercept, na.rm = T))

save(miss_delay_dur_LMM_FE, miss_delay_dur_LMM_RE, file = paste0(out_path,"LMM_state_data.RData"))
print("LMM finished")

#### Miss opp All ---------------------------------------------------------

get_miss_opp_res_GLMM <- function(trial_val){
  
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, inpatient, age_cat, female, rural, source, weekend, year, month,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, 
           #top2_high_inc_state_baddley
           state_abb) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  fit <- lme4::glmer(miss ~ inpatient + age_cat + female + rural + source + weekend + year + month + 
               asthma_prior_cp + copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp +(1|state_abb),
             data = reg_data,
             family = binomial(link = "logit"))
  
  RE_temp <- lme4::ranef(fit)$state_abb
  FE_temp <- lme4::fixef(fit)
  
  RE <- tibble(intercept = RE_temp[,1],
               state_abb = rownames(RE_temp))
  FE <- tibble(term = names(FE_temp),
               estimate = FE_temp)
  
  return(list(RE = RE,
              FE = FE))
  
}


miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                         function(x){get_miss_opp_res_GLMM(x)}, 
                                         mc.cores = num_cores)

RE <- do.call("rbind", lapply(miss_opp_res, function(x){x$RE} ))
FE <- do.call("rbind", lapply(miss_opp_res, function(x){x$FE} ))


miss_opp_GLMM_FE <- FE %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T) )

miss_opp_GLMM_RE <- RE %>% 
  group_by(state_abb) %>% 
  summarise(intercept_est = median(intercept, na.rm = T))

save(miss_opp_GLMM_FE, miss_opp_GLMM_RE, file = paste0(out_path,"GLMM_state_data.RData"))
print("GLMM finished")

################################################################################

library(tidyverse)
library(usmap)
out_path <- "/Volumes/Statepi_Diagnosis/projects/blasto/risk_models/"
load(paste0(out_path,"LMM_state_data.RData"))

states_w_info_LMM <- usmap::statepop %>% select(-pop_2015) %>%
  inner_join(miss_delay_dur_LMM_RE %>% rename(abbr = state_abb))

plot_usmap(data = states_w_info_LMM, values = "intercept_est")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name =  "",
                        label = scales::comma,) +
  labs(title = "LMM for Duration of Delay", subtitle = "Predicted State-level Random Intercept") +
  theme(legend.position = "right")



load(paste0(out_path,"GLMM_state_data.RData"))

states_w_info_GLMM <- usmap::statepop %>% select(-pop_2015) %>%
  inner_join(miss_opp_GLMM_RE %>% rename(abbr = state_abb))

plot_usmap(data = states_w_info_GLMM, values = "intercept_est")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name =  "",
                       label = scales::comma,) +
  labs(title = "GLMM for Log-odds of Miss Opportunities", subtitle = "Predicted State-level Random Intercept") +
  theme(legend.position = "right")


