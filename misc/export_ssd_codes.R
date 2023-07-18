
library(tidyverse)
library(codeBuildr)

codeBuildr::load_ssd_codes("sepsis_revised10") %>% 
  filter(!is.na(code)) %>% 
  rename(dx=code) %>% 
  select(dx) %>% 
  mutate(dx_ver=10) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sepsis/ssd_list/ssd_list.csv")


codeBuildr::load_ssd_codes("cvst") %>% 
  filter(!is.na(code)) %>% 
  rename(dx=code) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  select(dx,dx_ver,desc) %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/cvst/ssd_list/ssd_list.csv")

codeBuildr::load_ssd_codes("sarcoid") %>% 
  filter(!is.na(code)) %>% 
  rename(dx=code) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  select(dx,dx_ver,desc) %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/ssd_list/ssd_list.csv")
