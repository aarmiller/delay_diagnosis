
library(tidyverse)

load("/Shared/AML/truven_extracts/rx/all_abx/all_abx_dates.RData")

condition <- "histo"

# rx_dates %>% summarise(max(svcdate))

rx_dates <- rx_dates %>% distinct(enrolid,svcdate)

load(paste0("/Shared/AML/truven_extracts/dx/",condition,"/",condition,"_index_dx_dates.RData"))

index_dates <- index_dx_date %>% 
  filter(time_before_index>=365) %>% 
  filter(index_date<=18992) %>% 
  select(enrolid,index_date)

rx_counts <- index_dates %>% 
  left_join(rx_dates) %>% 
  mutate(days_since_index = svcdate-index_date) %>% 
  filter(between(days_since_index,-365,0)) %>% 
  distinct(enrolid,days_since_index) %>% 
  count(days_since_index)

save(rx_counts, file = paste0("/Shared/AML/tmp_transfer/abx_delay_grant/abx_before_data/",condition,".RData"))

rm(index_dates,rx_counts,condition)
gc()



condition_list <- c("blasto","cocci","pertussis","sepsis","tb","endocarditis",
                    "append","epidural_abs","hsv_enceph","nontb_myco","ebv","cmv",
                    "pe","chf","ami","sarcoid","dvt","interstitial_lung","venous_insuf",
                    "lung_cancer","bladder_cancer","non_hodgkins_lymphoma","hodgkins_lymphoma",
                    "renal_cell_cancer","acute_myeloid_leukemia","hairy_cell_leukemia",
                    "ulcerative_colitis","crohns","kawasaki","sle","wegeners","takayasu",
                    "periodic_fever","adult_stills","thyroiditis","behcet","ra","giant_cell_arteritis")

for (i in condition_list){
  
  print(i)
  
  load(paste0("/Shared/AML/truven_extracts/dx/",i,"/",i,"_index_dx_dates.RData"))
  
  index_dates <- index_dx_date %>% 
    filter(time_before_index>=365) %>% 
    filter(index_date<=18992) %>% 
    select(enrolid,index_date)
  
  rx_counts <- index_dates %>% 
    left_join(rx_dates) %>% 
    mutate(days_since_index = svcdate-index_date) %>% 
    filter(between(days_since_index,-365,0)) %>% 
    distinct(enrolid,days_since_index) %>% 
    count(days_since_index)
  
  save(rx_counts, file = paste0("/Shared/AML/tmp_transfer/abx_delay_grant/abx_before_data/",i,".RData"))
  
  rm(index_dates,rx_counts)
  gc()
  
}

codeBuildr::avail_disease_codes(F)

rm(list = ls())


condition_list <- c("histo","blasto","cocci","pertussis","sepsis","tb","endocarditis",
                    "append","epidural_abs","hsv_enceph","nontb_myco","ebv","cmv",
                    "pe","chf","ami","sarcoid","dvt","interstitial_lung","venous_insuf",
                    "lung_cancer","bladder_cancer","non_hodgkins_lymphoma","hodgkins_lymphoma",
                    "renal_cell_cancer","acute_myeloid_leukemia","hairy_cell_leukemia",
                    "ulcerative_colitis","crohns","kawasaki","sle","wegeners","takayasu",
                    "periodic_fever","adult_stills","thyroiditis","behcet","ra","giant_cell_arteritis")


plot_data <- tibble()

for (i in condition_list){
  
  load(paste0("/Volumes/AML/tmp_transfer/abx_delay_grant/abx_before_data/",i,".RData"))
  
  tmp <- rx_counts %>% 
    mutate(disease = i)
  
  plot_data <- bind_rows(plot_data,tmp)
  
  
}


for (i in condition_list){
  plot_data %>% 
    filter(disease==i) %>% 
    filter(days_since_index<0) %>% 
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    geom_smooth() +
    ggtitle(i)
  ggsave(paste0("~/Desktop/abx_plots/",i,".pdf"))
}

plot_data %>% 
  filter(disease==i) %>% 
  filter(days_since_index<0) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  geom_smooth() +
  ggtitle(i)
ggsave(paste0("~/Desktop/abx_plots/",i,".pdf"))
