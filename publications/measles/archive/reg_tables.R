

library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "measles"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 14 + 1
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
                          "REF", '  Age < 2', 
                         'age_cat(1,4]', '  Age 2-4'  ,
                         'age_cat(4,11]', '  Age 5-11'  , 
                         'age_cat(11,17]', '  Age 12-17'  , 
                         'age_cat(17,35]', '  Age 18-35'  , 
                         'age_cat(35,45]', '  Age 36-45'  , 
                         'age_cat(45,55]', '  Age 46-55'  ,
                         'age_cat(55,65]', '  Age 56-65'  , 
                         'age_cat(65,130]', '  Age >65'  , 
                         'femaleTRUE', 'Female'  , 
                         'rural', 'Rural location'  , 
                         "weekendTRUE", "Weekend Visit", 
                         'resp_antibiotic_drugs_window', 'Respiratory abx rx during delay window',
                         'vaccination_prior_cp', 'Measles vaccination prior to delay window'
                         )


month <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Month", "  January")),
                   tibble(term = paste0("month", 2:12),
                          label = paste0("  ", lubridate::month(2:12, label=TRUE, abbr = F))))

year <- bind_rows(tibble(term = c("Header", "REF"),
                          label = c("Year", "  /n2003")),
                   tibble(term = paste0("year", 2004:2023),
                          label = paste0("  /n", as.character(2004:2023))))

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
