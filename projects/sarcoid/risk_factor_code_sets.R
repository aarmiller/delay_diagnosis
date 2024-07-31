

## This script creates the risk factor code sets for the Sarcoid delay risk models.

library(tidyverse)

###################
#### Load Data ####
###################


### Load Redbook ---------------------------------------------------------------

redbook <- read_csv("/Shared/Statepi_Marketscan/databases/Truven/medication_dbs/redbook.csv")
# redbook <- read_csv("/Volumes/Statepi_Marketscan/databases/Truven/medication_dbs/redbook.csv")


### Load Codes from Miles ------------------------------------------------------

tmp_inhaler <- readxl::read_xlsx("/Shared/Statepi_Diagnosis/projects/sarcoid/external_code_sets/inhaler_codes MH.xlsx")
tmp_steroid <- readxl::read_xlsx("/Shared/Statepi_Diagnosis/projects/sarcoid/external_code_sets/oral_steroid_codes MH.xlsx")
tmp_ther_class <- readxl::read_xlsx("/Shared/Statepi_Diagnosis/projects/sarcoid/external_code_sets/redbook_therapeutic_classes MH.xlsx")
tmp_ther_group <- readxl::read_xlsx("/Shared/Statepi_Diagnosis/projects/sarcoid/external_code_sets/redbook_therapeutic_detailed_group MH.xlsx")

names(tmp_inhaler) <- tolower(names(tmp_inhaler))
names(tmp_steroid) <- tolower(names(tmp_steroid))
names(tmp_ther_class) <- tolower(names(tmp_ther_class))
names(tmp_ther_group) <- tolower(names(tmp_ther_group))

tmp_inhaler <- tmp_inhaler %>% 
  mutate(ndcnum = str_pad(ndcnum,width = 11,side = "left", pad = "0"))

tmp_steroid <- tmp_steroid %>% 
  mutate(ndcnum = str_pad(ndcnum,width = 11,side = "left", pad = "0"))


############################
#### Export Final Codes ####
############################

### Oral Steroids --------------------------------------------------------------

steroids <- tmp_steroid %>% 
  filter(maintds %in% c("Both acute and chronic","Primarily acute")) 

steroids <- steroids %>% 
  distinct(gennme) %>% 
  inner_join(redbook) %>% 
  filter(roads=="Oral") %>% 
  distinct(ndcnum)


#### Inhalers ------------------------------------------------------------------

inhalers <- tmp_inhaler %>%
  distinct(gennme) %>% 
  inner_join(redbook) %>% 
  filter(roads=="Inhalation") %>% 
  distinct(ndcnum)


### new codes to add from therapuetic group ------------------------------------

add_codes <- tmp_ther_group %>% 
  select(thrdtds,drug_class=`drug class`) %>% 
  filter(!is.na(drug_class)) %>% 
  mutate(drug_class = tolower(drug_class)) %>% 
  mutate(drug_class = ifelse(drug_class=="antibiotic (tb related)","antibiotic",drug_class)) %>% 
  inner_join(redbook) %>% 
  select(drug_class,ndcnum)

tmp_inhaler_add <- add_codes %>% 
  filter(drug_class=="inhaler") %>% 
  inner_join(redbook) %>% 
  filter(roads == "Inhalation") %>% 
  distinct(ndcnum)

tmp_steroid_add <- add_codes %>% 
  filter(drug_class=="steroid") %>% 
  inner_join(redbook) %>% 
  filter(roads == "Oral") %>% 
  distinct(ndcnum)

## finalize --------------------------------------------------------------------


steroids <- bind_rows(steroids,tmp_steroid_add) %>% distinct()
inhalers <- bind_rows(inhalers,tmp_inhaler_add) %>% distinct()

other_drugs <- add_codes %>% 
  filter(!(drug_class %in% c("inhaler","steroid")))




save(steroids,inhalers,other_drugs, file = "/Shared/Statepi_Diagnosis/projects/sarcoid/code_sets/drug_risk_factor_codes.RData")



