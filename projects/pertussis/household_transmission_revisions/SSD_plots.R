

rm(list = ls())

library(tidyverse)


db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pertussis/pertussis.db")

tmp <- db %>% 
  tbl("all_proc_visits") %>% 
  filter(proc %in% c("87798","86615")) %>% 
  filter(between(days_since_index,-14,14)) %>% 
  collect()

# Validated ids
validated_ids <- tmp %>% 
  distinct(patient_id)
rm(tmp)

# Exclude Medicaid
include_ids <- db %>% 
  tbl("enrolid_crosswalk") %>% 
  filter(medicaid==0) %>% 
  select(patient_id) %>% 
  collect()

dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(include_ids, copy = TRUE) %>% 
  collect()

dx_visits <- dx_visits %>% 
  left_join(mutate(validated_ids,validated=1L)) %>% 
  mutate(validated = replace_na(validated,0L))

save(dx_visits, file = "/Shared/AML/tmp_transfer/pertussis_dx_visits.RData")

#### RUN LOCALLY ####
rm(list = ls())


load("/Volumes/AML/tmp_transfer/pertussis_dx_visits.RData")

# Load New SSD set
ssd_codes <- read_csv("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/revised_ssd.csv")

ssd_codes <- ssd_codes %>% 
  filter(Include!="No") %>% 
  select(dx,dx_ver,desc)

ssd_counts <- dx_visits %>% 
  filter(between(days_since_index,-180,180)) %>% 
  inner_join(ssd_codes) %>% 
  distinct(patient_id,dx,dx_ver,days_since_index) %>% 
  count(dx,dx_ver,days_since_index)

ssd_counts_validated <- dx_visits %>% 
  filter(validated==1) %>% 
  filter(between(days_since_index,-180,180)) %>% 
  inner_join(ssd_codes) %>% 
  distinct(patient_id,dx,dx_ver,days_since_index) %>% 
  count(dx,dx_ver,days_since_index)


ssd_counts <- ssd_counts %>% 
  inner_join(codeBuildr::all_icd_labels)

ssd_counts_validated <- ssd_counts_validated %>% 
  inner_join(codeBuildr::all_icd_labels)

ssd_counts_totals <- ssd_counts %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

ssd_counts_totals_validated <- ssd_counts_validated %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

ssd_counts <- ssd_counts %>% 
  mutate(label = paste0("ICD-",dx_ver," ",dx,": ",str_sub(desc,1,40)))

ssd_counts_validated <- ssd_counts_validated %>% 
  mutate(label = paste0("ICD-",dx_ver," ",dx,": ",desc))

ssd_counts_totals

max(ssd_counts_totals$group)

ssd_counts_totals <- ssd_counts_totals %>% 
  mutate(row = row_number(),
         group = (row-7)%/%8+2)

pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/ssd_trends.pdf",onefile = TRUE)
pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/final_revisions/sdc2.pdf",onefile = TRUE)
for (i in 1:max(ssd_counts_totals$group)){
  
  if (i==1){
    p <- ssd_counts_totals %>% 
      filter(group==i) %>%  
      select(dx,dx_ver) %>% 
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
      theme_bw() +
      ggtitle("Supplemental Digital Content 2") +
      ylab("Number of Visits with Diagnosis") +
      xlab("Weeks Since Index")
  } else {
    p <- ssd_counts_totals %>% 
      filter(group==i) %>%  
      select(dx,dx_ver) %>% 
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
      theme_bw() +
      ylab("Number of Visits with Diagnosis") +
      xlab("Weeks Since Index")
  }
  
  print(p)
}
dev.off()

for (i in 2:max(ssd_counts_totals$group)){
  pdf(paste0("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/final_revisions/sdc2/set_",i,".pdf"),height = 9.5)
  p <- ssd_counts_totals %>% 
    filter(group==i) %>%  
    select(dx,dx_ver) %>% 
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
    theme_bw() +
    ylab("Number of Visits with Diagnosis") +
    xlab("Weeks Since Index")
  print(p)
  dev.off()
}

i <- 11
pdf(paste0("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/final_revisions/sdc2/set_",i,".pdf"),height = 3)
p <- ssd_counts_totals %>% 
  filter(group==i) %>%  
  select(dx,dx_ver) %>% 
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
  theme_bw() +
  ylab("Number of Visits with Diagnosis") +
  xlab("Weeks Since Index")
print(p)
dev.off()

 # ssd_counts_totals %>% 
#   slice((i*6-5):(i*6))
# 
# ssd_counts_totals %>% 
#   inner_join(distinct(ssd_counts,dx,dx_ver,desc)) %>% 
#   write_csv("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/ssd_include.csv")

pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/ssd_trends_validated.pdf",onefile = TRUE)
for (i in 1:ceiling(nrow(ssd_counts_totals_validated)/6)){
  
  p <- ssd_counts_totals_validated %>% 
    slice((i*6-5):(i*6)) %>% 
    select(-n) %>% 
    inner_join(ssd_counts_validated) %>% 
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

#######################
#### Plot Selected ####
#######################

count_weeks <- function(codes){
  dx_visits %>% 
    filter(dx %in% c(codes)) %>% 
    filter(between(days_since_index,-180,180)) %>% 
    filter(days_since_index!=0) %>% 
    mutate(weeks_since_index = days_since_index %/% 7) %>% 
    distinct(patient_id,weeks_since_index) %>% 
    count(weeks_since_index)
}

count_weeks2 <- function(codes){
  dx_visits %>% 
    filter(dx %in% c(codes)) %>% 
    filter(between(days_since_index,-180,180)) %>% 
    filter(days_since_index!=0) %>% 
    distinct(patient_id,days_since_index) %>% 
    count(days_since_index) %>% 
    mutate(weeks_since_index = days_since_index %/% 7) %>% 
    group_by(weeks_since_index) %>% 
    summarise(n = mean(n,na.rm = T))
}

# cough
count_weeks(c("7862","R05")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

# fever
count_weeks(c("7806","R509")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

# croup
count_weeks(c("4644","J050")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

# Asthma
count_weeks(c("49392","J45901")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

# Bronchitis
count_weeks2(c("4660","J209")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

# Unspecified resp infection
count_weeks2(c("4659","J069")) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() 

pdf("~/OneDrive - University of Iowa/WorkingPapers/alan_pertussis/delays and household transmisson/submissions/PIDJ/revisions/focal_ssd_trends.pdf",
    width = 6, height = 4)
bind_rows(mutate(count_weeks2(c("4659","J069")),
                 group = "Acute Resp. Inf. Unspec."),
          mutate(count_weeks2(c("4660","J209")),
                 group = "Bronchitis"),
          mutate(count_weeks2(c("49392","J45901")),
                 group = "Asthma"),
          mutate(count_weeks2(c("4644","J050")),
                 group = "Croup"),
          mutate(count_weeks2(c("7806","R509")),
                 group = "Fever"),
          mutate(count_weeks2(c("7862","R05")),
                 group = "Cough"))   %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal() +
  facet_wrap(~group,scales = "free_y") +
  xlab("Weeks Since Index Pertussis Diagnosis") +
  ylab("Number of Patients with a Visit")
dev.off()


dx_visits %>% 
  filter(dx %in% c("7862","R05")) %>% 
  filter(between(days_since_index,-180,180)) %>% 
  filter(days_since_index!=0) %>% 
  mutate(weeks_since_index = days_since_index %/% 7) %>% 
  distinct(patient_id,weeks_since_index) %>% 
  count(weeks_since_index) %>% 
  ggplot(aes(weeks_since_index,n)) +
  geom_point() +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal()



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
