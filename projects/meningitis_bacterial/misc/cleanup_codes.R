


mening_codes_all <- codeBuildr::load_disease_codes("meningitis")$meningitis

prior_list <- read_csv("~/Documents/GitHub/delay_dx/params/index_condition_codes/index_meningitis.csv")

prior_list %>% anti_join(bind_rows(tibble(code = mening_codes_all$icd9_codes,icd_version = 9),
                                    tibble(code = mening_codes_all$icd10_codes,icd_version = 10)))

tmp <- bind_rows(tibble(dx = mening_codes_all$icd9_codes,dx_ver = 9),
          tibble(dx = mening_codes_all$icd10_codes,dx_ver = 10)) %>% left_join(codeBuildr::all_icd_labels)

tmp %>% filter(!str_detect(tolower(desc),"meningitis"))

prior_list %>% rename(dx=code,dx_ver=icd_version) %>% left_join(codeBuildr::all_icd_labels)

codeBuildr::all_icd_labels %>% filter(dx=="360")

bind_rows(tibble(dx = mening_codes_all$icd9_codes,icd_version = 9),
          tibble(dx = mening_codes_all$icd10_codes,icd_version = 10))

codeBuildr::all_icd_labels %>% filter(dx_ver==9) %>% 
  inner_join(tibble(dx = tmp$meningitis$icd9_codes)) %>% 
  distinct()


tmp$meningitis$icd9_codes

library(tidyverse)
db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/meningitis/meningitis.db")

db %>% tbl("all_rx_visits")

rm(db)

# Cleanup validation procedure codes
library(tidyverse)

tmp <- read_csv("/Volumes/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/top_1000_labels_cleaned.csv")

tmp %>% filter(is.na(code_type)) %>% write_csv("/Volumes/Statepi_Diagnosis/projects/meningitis_bacterial/prelim_data/codes_to_label.csv")
