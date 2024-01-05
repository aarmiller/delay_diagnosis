

rm(list = ls())

library(tidyverse)


db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis.db")

tmp <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("87798","86615")) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect()

include_ids <- tmp %>% 
  distinct(patient_id)

dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(include_ids, copy = TRUE) %>% 
  collect()

save(dx_visits, file = "/Shared/AML/tmp_transfer/pertussis_dx_visits.RData")

load("/Volumes/AML/tmp_transfer/pertussis_dx_visits.RData")

# Load SSD set
ssd_codes <- read_csv("~/Documents/GitHub/delay_dx/params/ssd_codes/pertus/ssd_codes.csv")

ssd_codes <- ssd_codes %>% 
  select(dx = icd_codes, dx_ver = icd_version)

ssd_counts <- dx_visits %>% 
  filter(between(days_since_index,-180,180)) %>% 
  inner_join(ssd_codes) %>% 
  distinct(patient_id,dx,dx_ver,days_since_index) %>% 
  count(dx,dx_ver,days_since_index)


ssd_counts %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n))

ssd_counts <- ssd_counts %>% 
  inner_join(codeBuildr::all_icd_labels)

ssd_counts_totals <- ssd_counts %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

ssd_counts <- ssd_counts %>% 
  mutate(label = paste0("ICD-",dx_ver," ",dx,": ",desc))


pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/ssd_trends.pdf",onefile = TRUE)
for (i in 1:ceiling(nrow(ssd_counts_totals)/6)){
  
  p <- ssd_counts_totals %>% 
    slice((i*6-5):(i*6)) %>% 
    select(-n) %>% 
    inner_join(ssd_counts) %>% 
    filter(days_since_index!=0) %>% 
    group_by(dx,dx_ver,label) %>% 
    mutate(weeks_since_index = days_since_index %/% 7) %>% 
    group_by(dx,dx_ver,label,weeks_since_index) %>% 
    summarise(n = mean(n)) %>% 
    ggplot(aes(weeks_since_index,n)) +
    geom_point() +
    facet_wrap(~label, scales = "free_y", ncol = 2) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    theme_bw()
  
  print(p)
}
dev.off()

ssd_counts_totals %>% 
  slice((i*6-5):(i*6))

ssd_counts_totals %>% 
  inner_join(distinct(ssd_counts,dx,dx_ver,desc)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/ssd_include.csv")
  


ssd_counts_totals %>% 
  slice(1:6) %>% 
  select(-n) %>% 
  inner_join(ssd_counts) %>% 
  filter(days_since_index!=0) %>% 
  group_by(dx,dx_ver,label) %>% 
  mutate(weeks_since_index = days_since_index %/% 7) %>% 
  group_by(dx,dx_ver,label,weeks_since_index) %>% 
  summarise(n = mean(n)) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  facet_wrap(~label, scales = "free_y", ncol = 2) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_bw()



### Analyze Cluster Results ###
rm(list = ls())
source("cluster/cluster_functions.R")

load("/Volumes/Statepi_Diagnosis/prelim_results/pertussis/cluster_results/cluster_res.RData")

plot_cluster_shapes(clust_res_10,10)

focal9 <- find_focal_clusters(clust_res_9,focal_code = "7862",dx_ver = 9)

focal10 <- find_focal_clusters(clust_res_10,focal_code = "R05",dx_ver = 10)

ssd_codes <- read_csv("~/Documents/GitHub/delay_dx/params/ssd_codes/pertus/ssd_codes.csv")

ssd_codes <- ssd_codes %>% 
  select(dx = icd_codes, dx_ver = icd_version)

ssd_codes %>% 
  filter(dx_ver==9) %>% 
  mutate(included_ssd=1L) %>% 
  select(-dx_ver)

focal9_include <- focal9$condition_counts %>% 
  filter(n>1) %>% 
  left_join(ssd_codes %>% 
             filter(dx_ver==9) %>% 
             select(code= dx,-dx_ver) %>% 
              mutate(included_ssd=1L)) %>% 
  mutate(included_ssd=replace_na(included_ssd,0L)) %>% 
  distinct() %>% 
  mutate(page_no = ((row_number()-1) %/% 6)+1) %>% 
  mutate(index = row_number()) %>% 
  select(index,page_no,everything())

focal9_include %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/cluster_res_icd9.csv")


focal10_include <- focal10$condition_counts %>% 
  filter(n>1) %>% 
  left_join(ssd_codes %>% 
              filter(dx_ver==10) %>% 
              select(code= dx,-dx_ver) %>% 
              mutate(included_ssd=1L)) %>% 
  mutate(included_ssd=replace_na(included_ssd,0L)) %>% 
  distinct() %>% 
  mutate(page_no = ((row_number()-1) %/% 6)+1) %>% 
  mutate(index = row_number()) %>% 
  select(index,page_no,everything())
  
focal10_include %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/cluster_res_icd10.csv")



load("/Volumes/AML/tmp_transfer/pertussis_dx_visits.RData")

select(focal9_include,index,dx=code)

tmp1 <- dx_visits %>% 
  filter(between(days_since_index,-180,180)) %>% 
  filter(dx_ver==9) %>% 
  inner_join(select(focal9_include,index,dx=code)) %>% 
  distinct(patient_id,dx,index,days_since_index) %>% 
  count(dx,index,days_since_index) 



pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/clust_plots_icd9.pdf",onefile = TRUE)
for (i in 1:max(focal9_include$page_no)){
  
  p <- focal9_include %>% 
    filter(page_no==i) %>% 
    select(dx=code,desc) %>% 
    inner_join(tmp1) %>% 
    filter(days_since_index!=0) %>% 
    mutate(label = paste0(index," - ICD-9", dx, ": ",str_sub(desc,1,30))) %>% 
    group_by(dx,label) %>% 
    mutate(weeks_since_index = days_since_index %/% 7) %>% 
    group_by(dx,label,weeks_since_index) %>% 
    summarise(n = mean(n)) %>% 
    ggplot(aes(weeks_since_index,n)) +
    geom_point() +
    facet_wrap(~label, scales = "free_y", ncol = 2) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    theme_bw()
  
  print(p)
}
dev.off()



tmp1 <- dx_visits %>% 
  filter(between(days_since_index,-180,180)) %>% 
  filter(dx_ver==10) %>% 
  inner_join(select(focal10_include,index,dx=code)) %>% 
  distinct(patient_id,dx,index,days_since_index) %>% 
  count(dx,index,days_since_index) 



pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/clust_plots_icd10.pdf",onefile = TRUE)
for (i in 1:max(focal10_include$page_no)){
  
  p <- focal10_include %>% 
    filter(page_no==i) %>% 
    select(dx=code,desc) %>% 
    inner_join(tmp1) %>% 
    filter(days_since_index!=0) %>% 
    mutate(label = paste0(index," - ICD-10", dx, ": ",str_sub(desc,1,30))) %>% 
    group_by(dx,label) %>% 
    mutate(weeks_since_index = days_since_index %/% 7) %>% 
    group_by(dx,label,weeks_since_index) %>% 
    summarise(n = mean(n)) %>% 
    ggplot(aes(weeks_since_index,n)) +
    geom_point() +
    facet_wrap(~label, scales = "free_y", ncol = 2) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    theme_bw()
  
  print(p)
}
dev.off()
