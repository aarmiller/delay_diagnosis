

library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "sarcoid"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

load(paste0(out_path,"ssd_miss_risk_models.RData"))

table_out_path <- paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/reg_tables/")

# setup better labels
better_labels <- tribble(~term,~label, 
                          "Header", 'Database Source',
                          "REF", '  Commercial', 
                         'sourcemdcr', '  Medicare'  , 
                         'sourcemedicaid', '  Medicaid', 
                          "Header", 'Age', 
                          "REF", '  Age < 18', 
                         'age_cat(17,34]', '  Age 18-34'  ,
                         'age_cat(34,44]', '  Age 35-44'  , 
                         'age_cat(44,54]', '  Age 45-54'  , 
                         'age_cat(54,64]', '  Age 55-64'  , 
                         'age_cat(64,130]', '  Age >= 65'  , 
                         'femaleTRUE', 'Female'  , 
                         "any_obesity", "Any obesity dx prior to index", 
                         "inpatient", "Inpatient", 
                         'rural', 'Rural location'  , 
                         # 'setting_labelAll three', 'Outpatient, ED, and Inpatient'  ,
                         # 'setting_labelED and inpatient', 'ED and Inpatient'  ,
                         # 'setting_labelED only', 'ED Only'  ,
                         # 'setting_labelInpatient only', 'Inpatient Only'  ,
                         # 'setting_labelOut and ED', 'Outpatient and ED'  ,
                         # 'setting_labelOut and inpatient', 'Outpatient and Inpatient'  ,
                         # 'setting_labelnone', 'Not Outpatient, ED or Inpatient'  ,
                         "weekendTRUE", "Weekend Visit", 
                         "antiacid_ppi_drugs_window", "Antacid/PPI rx during delay window", 
                         "antibiotic_drugs_window", "Antibiotic rx during delay window", 
                         "antihistamine_drugs_window", "Antihistamine rx during delay window", 
                         "cough_suppressant_drugs_window", "Cough suppressant rx during delay window", 
                         "diuretic_drugs_window", "Diuretic rx during delay window", 
                         "inhalers_window", "Inhaler rx during delay window", 
                         "nasal_spray_drugs_window", "Nasal spray rx during delay window", 
                         "oral_steroids_window", "Oral steroid rx during delay window")


month <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Month", "  January")),
                   tibble(term = paste0("month", 2:12),
                          label = paste0("  ", lubridate::month(2:12, label=TRUE, abbr = F))))

year <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Year", "  /n2003")),
                   tibble(term = paste0("year", 2004:2022),
                          label = paste0("  /n", as.character(2004:2022))))

master_table <- bind_rows(better_labels, month, year)

#############
## Table 1 ##
#############

data <- ssd_miss_risk_models$miss_opp_res_inpatient_ind

miss_opp_table <- left_join(master_table, data) %>% 
  filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
  mutate("OR (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
                                               " (",
                                               trimws(format(round(low, 2), nsmall = 2)),
                                               "-",
                                               trimws(format(round(high, 2), nsmall = 2)),
                                               ")"),
                                               ifelse(term == "Header",  "", term))) %>% 
  select(2, 6)

write_csv(x = miss_opp_table, file = paste0(table_out_path,"miss_opp_table.csv") )


#############
## Table 2 ##
#############

data <- ssd_miss_risk_models$miss_delay_pat_res

miss_delay_pat_res <- left_join(master_table, data) %>% 
  filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
  mutate("OR (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
                                                   " (",
                                                   trimws(format(round(low, 2), nsmall = 2)),
                                                   "-",
                                                   trimws(format(round(high, 2), nsmall = 2)),
                                                   ")"),
                               ifelse(term == "Header",  "", term))) %>% 
  select(2, 6)

write_csv(x = miss_delay_pat_res, file = paste0(table_out_path,"miss_delay_pat_res.csv") )


#############
## Table 3 ##
#############

data <- ssd_miss_risk_models$miss_dur_res_weibull

miss_dur_res <- left_join(master_table, data) %>% 
  filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
  mutate("Estimate (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
                                                   " (",
                                                   trimws(format(round(low, 2), nsmall = 2)),
                                                   "-",
                                                   trimws(format(round(high, 2), nsmall = 2)),
                                                   ")"),
                               ifelse(term == "Header",  "", term))) %>% 
  select(2, 6)

write_csv(x = miss_dur_res, file = paste0(table_out_path,"miss_dur_res.csv") )

#############
## Table 4 ##
#############

treatment_tbl <- tribble(~term,~label, 
                         "antiacid_ppi_drugs_window", "Antacid/PPI rx during delay window", 
                         "antibiotic_drugs_window", "Antibiotic rx during delay window", 
                         "antihistamine_drugs_window", "Antihistamine rx during delay window", 
                         "cough_suppressant_drugs_window", "Cough suppressant rx during delay window", 
                         "diuretic_drugs_window", "Diuretic rx during delay window", 
                         "inhalers_window", "Inhaler rx during delay window", 
                         "nasal_spray_drugs_window", "Nasal spray rx during delay window", 
                         "oral_steroids_window", "Oral steroid rx during delay window")

treatment_res <- miss_delay_pat_res %>% inner_join(treatment_tbl) %>% 
  inner_join(miss_dur_res) %>% 
  select(-term)

names(treatment_res) <- c("Treatment", "Delayed Patient Model\n OR (95%CI)", "Duration of Delay Model\n Estimate (95%CI)")
write_csv(x = treatment_res, file = paste0(table_out_path,"treatment_res.csv") )
