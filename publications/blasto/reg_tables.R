

library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "blasto"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 63 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <- paste0(delay_params$out_path,"risk_models/") 

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
                         'rural', 'Rural location'  , 
                         "weekendTRUE", "Weekend Visit", 
                         "asthma_prior_cp", "Asthma dx prior to delay window", 
                         "copd_prior_cp", "COPD dx prior to delay window", 
                         "chest_ct_prior_cp", "Chest CT prior to delay window", 
                         "chest_xray_prior_cp", "Chest X-ray prior to delay window", 
                         "Est ID_consult vs no ID_consult - top2_high_inc_state_baddley", "Infectious Disease consult during visit - Residing in top 2 highest blastomycosis incidence state",
                         "Est ID_consult vs no ID_consult - Non-top2_high_inc_state_baddley",  "Infectious Disease consult during visit - Not residing in top 2 highest blastomycosis incidence state"
                         )


month <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Month", "  January")),
                   tibble(term = paste0("month", 2:12),
                          label = paste0("  ", lubridate::month(2:12, label=TRUE, abbr = F))))

year <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Year", "  /n2002")),
                   tibble(term = paste0("year", 2003:2022),
                          label = paste0("  /n", as.character(2003:2022))))

master_table <- bind_rows(better_labels, month, year)

#############
## Table 1 ##
#############

data <- ssd_miss_risk_models$miss_opp_res_interaction

miss_opp_table <- left_join(master_table, data) %>% 
  filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
  mutate("OR (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
                                               " (",
                                               trimws(format(round(low, 2), nsmall = 2)),
                                               "-",
                                               trimws(format(round(high, 2), nsmall = 2)),
                                               ")"),
                                               ifelse(term == "Header",  "", term))) %>% 
  select(2, 8)

write_csv(x = miss_opp_table, file = paste0(table_out_path,"miss_opp_table.csv") )


#############
## Table 2 ##
#############

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
                         'rural', 'Rural location'  , 
                         "weekendTRUE", "Weekend Visit", 
                         "asthma_prior_cp", "Asthma dx prior to delay window", 
                         "copd_prior_cp", "COPD dx prior to delay window", 
                         "chest_ct_prior_cp", "Chest CT prior to delay window", 
                         "chest_xray_prior_cp", "Chest X-ray prior to delay window", 
                         'resp_antibiotic_drugs_window', 'Respiratory abx rx during delay window',
                         'inhalers_window', 'Inhaler rx during delay window',
                         'top2_high_inc_state_baddley', 'Residing in top 2 highest blastomycosis incidence state',
)
master_table <- bind_rows(better_labels, month, year)

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
  select(2, 8)

write_csv(x = miss_delay_pat_res, file = paste0(table_out_path,"miss_delay_pat_res.csv") )


#############
## Table 3 ##
#############

data <- ssd_miss_risk_models$miss_dur_res

miss_dur_res <- left_join(master_table, data) %>% 
  filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
  mutate("Estimate (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
                                                   " (",
                                                   trimws(format(round(low, 2), nsmall = 2)),
                                                   "-",
                                                   trimws(format(round(high, 2), nsmall = 2)),
                                                   ")"),
                               ifelse(term == "Header",  "", term))) %>% 
  select(2, 8)

write_csv(x = miss_dur_res, file = paste0(table_out_path,"miss_dur_res.csv") )
