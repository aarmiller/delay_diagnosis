
library(tidyverse)
library(lubridate)
library(smallDB)

# devtools::install_github("aarmiller/smallDB")

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/blasto/table_fun.R")

## tables for full cohort

full_cohort <- main_fun("blasto", cp = 63)

## tables for inpatient cohort

inpatient_cohort <- main_fun("blasto_inpatient", cp = 63)

## tables for outpatient cohort

outpatient_cohort <- main_fun("blasto_outpatient", cp = 63)


# Table 1 ----------------------------------------------------------------------

full_cohort$table1 %>% 
  rename('Full Cohort' = 'Total Patients (% of patients)') %>% 
  inner_join(inpatient_cohort$table1 %>% 
               rename('Inpatient Cohort' = 'Total Patients (% of patients)'), by = "Characteristic") %>% 
  inner_join(outpatient_cohort$table1 %>% 
               rename('Outpatient Cohort' = 'Total Patients (% of patients)'), by = "Characteristic") %>% 
  write_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", "blasto", "/tables/table1_combined.csv"))

# Table 2 ----------------------------------------------------------------------

table2_combined <- full_cohort$table2 %>% 
  mutate(label = 'Full') %>% 
  bind_rows(inpatient_cohort$table2 %>% 
              mutate(label = 'Inpatient')) %>% 
  bind_rows(outpatient_cohort$table2 %>% 
              mutate(label = 'Outpatient')) %>% 
  mutate(Setting = factor(Setting, levels = c("outpatient", "ed", "obs_stay", "inpatient", "inpatient visit"),
                          labels = c("Oupatient", "ED", "Observational Stay", "Inpatient", "Inpatient visit"))) %>% 
  mutate(label = format(label, levels = c('Full', 'Inpatient', 'Outpatient')))

table2_combined %>% select(Setting, label, everything()) %>% 
  arrange(Setting, label) %>% 
  write_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", "blasto", "/tables/table2_combined.csv"))


# Table 3 ----------------------------------------------------------------------

full_cohort$table3 %>% 
  rename(estimate_Full_Cohort = estimate,
         CI_Full_Cohort = CI) %>% 
  inner_join(inpatient_cohort$table3 %>% 
              rename(estimate_Inpatient = estimate,
                     CI_Inpatient = CI)) %>% 
  inner_join(outpatient_cohort$table3 %>% 
               rename(estimate_Outpatient = estimate,
                      CI_Outpatient = CI)) %>% 
  write_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", "blasto", "/tables/table3_combined.csv"))


# # Table 4 ----------------------------------------------------------------------
# 
# better_labels <- tribble(~term,~label, 
#                          "Header", 'Database Source',
#                          "REF", '  Commercial', 
#                          'sourcemdcr', '  Medicare'  , 
#                          'sourcemedicaid', '  Medicaid', 
#                          "Header", 'Age', 
#                          "REF", '  Age < 18', 
#                          'age_cat(17,34]', '  Age 18-34'  ,
#                          'age_cat(34,44]', '  Age 35-44'  , 
#                          'age_cat(44,54]', '  Age 45-54'  , 
#                          'age_cat(54,64]', '  Age 55-64'  , 
#                          'age_cat(64,130]', '  Age >= 65'  , 
#                          'femaleTRUE', 'Female'  , 
#                          "weekendTRUE", "Weekend Visit",
#                          'abx', 'Antibiotic prescription during visit',
#                          "opioid", "Opioid prescription during visit", 
#                          "ID_consult", "Infectious Disease consult during visit"
# )
# 
# reg_full <- read_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/blasto/tables/table4.csv"))
# reg_validated <- read_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/blasto/tables/table4_validated.csv"))
# 
# left_join(better_labels, reg_full) %>% 
#   filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
#   mutate("FC OR (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
#                                                    " (",
#                                                    trimws(format(round(low, 2), nsmall = 2)),
#                                                    "-",
#                                                    trimws(format(round(high, 2), nsmall = 2)),
#                                                    ")"),
#                                ifelse(term == "Header",  "", term))) %>% 
#   select(2, 7) %>% inner_join(
#     left_join(better_labels, reg_validated) %>% 
#       filter(!is.na(est) | (term %in% c("Header", "REF")) ) %>% 
#       mutate("VC OR (95%CI)" = ifelse(!is.na(est), paste0(trimws(format(round(est, 2), nsmall = 2)), 
#                                                        " (",
#                                                        trimws(format(round(low, 2), nsmall = 2)),
#                                                        "-",
#                                                        trimws(format(round(high, 2), nsmall = 2)),
#                                                        ")"),
#                                    ifelse(term == "Header",  "", term))) %>% 
#       select(2, 7)
#   ) %>% 
#   write_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", "blasto", "/tables/table4_combined.csv"))



# Appendix table ----------------------------------------------------------------------

app_combined <- full_cohort$table4 %>% 
  mutate(label = 'Full') %>% 
  bind_rows(inpatient_cohort$table4 %>% 
              mutate(label = 'Inpatient')) %>% 
  bind_rows(outpatient_cohort$table4 %>% 
              mutate(label = 'Outpatient')) %>% 
  mutate(metric = factor(metric, levels = c("percent_patient_delayed", "mean_n_miss", "mean_dur"),
                          labels = c("Percent of patients that experienced a delay (95%CI)", "Mean number of missed opportunities per patient (95%CI)", "Mean delay duration (days) (95%CI)"))) %>% 
  mutate(label = format(label, levels = c('Full', 'Inpatient','Outpatient')))

app_combined %>% select(metric, label, everything()) %>% 
  arrange(metric, label) %>% 
  write_csv(paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", "blasto", "/tables/app_combined.csv"))
