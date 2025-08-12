rm(list = ls())

library(tidyverse)
library(bit64)

################
#### params ####
################

age_breaks <- c(-1,17,40,120)


###################
#### Load Data ####
###################


load("/Volumes/Statepi_Diagnosis/projects/cocci/index_cases.RData")

index_cases 


load("/Volumes/Statepi_Diagnosis/prelim_results/cocci/delay_results/all_dx_visits.RData")

ssd_codes <- codeBuildr::load_ssd_codes("cocci") %>% 
  mutate(dx_ver = ifelse(type == "icd9",9L,10L)) %>% 
  select(dx=code,dx_ver)

ssd_dx_visits <- all_dx_visits %>% 
  inner_join(select(index_cases,patient_id)) %>% 
  inner_join(ssd_codes)

load("/Volumes/Statepi_Diagnosis/prelim_results/cocci/delay_results/demo_data.RData")

demo <- demo1 %>% 
  mutate(age = year(as_date(index_date))-dobyr) %>% 
  select(patient_id,age,sex)

age_cat_counts <- demo %>% 
  mutate(age_cat = cut(age,breaks = age_breaks)) %>% 
  count(age_cat,name = "age_count") 

ssd_dx_visits <- ssd_dx_visits %>% 
  inner_join(demo) %>% 
  mutate(age_cat = cut(age,breaks = age_breaks)) %>% 
  mutate(week = -(days_since_index %/% 7)) 


dx_week_counts <- ssd_dx_visits %>% 
  filter(between(week,1,52)) %>% 
  count(dx,dx_ver,age_cat,week,dx) %>% 
  inner_join(age_cat_counts) %>% 
  mutate(frac = 100*n/age_count) %>% 
  inner_join(codeBuildr::all_icd_labels)

dx_week_counts <- dx_week_counts %>%
  mutate(label = paste0("dx",dx_ver," ",dx,": ",str_sub(desc,1,40)))

############################################
##### Counts in and out of delay window ####
############################################

tmp1 <- ssd_dx_visits %>% 
  filter(days_since_index<0) %>% 
  mutate(period = ifelse(days_since_index < -13*7,"before","during")) %>% 
  count(dx,dx_ver,age_cat,period) %>% 
  inner_join(age_cat_counts) %>% 
  mutate(frac = round(100*n/age_count,3)) %>% 
  select(dx:period,frac) 


tmp2 <- ssd_dx_visits %>% 
  filter(days_since_index<0) %>% 
  mutate(period = ifelse(days_since_index < -13*7,"before","during")) %>% 
  count(dx,dx_ver,period) %>% 
  mutate(frac = round(100*n/sum(age_cat_counts$age_count),3)) %>% 
  select(dx:period,frac) 

during_counts <- tmp1 %>% 
  filter(period == "during") %>% 
  select(-period) %>% 
  spread(key = age_cat, value = frac) %>% 
  mutate_all(~replace_na(.,0L)) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  select(dx,dx_ver,desc,"<18"=`(-1,17]`,"18-40"=`(17,40]`,">40"=`(40,120]`) %>% 
  left_join(tmp2 %>% 
              filter(period == "during") %>% 
              select(dx,dx_ver,"All Ages"=frac)) %>% 
  arrange(desc(`All Ages`))
  
before_counts <- tmp1 %>% 
  filter(period == "before") %>% 
  select(-period) %>% 
  spread(key = age_cat, value = frac) %>% 
  mutate_all(~replace_na(.,0L)) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  select(dx,dx_ver,desc,"before <18"=`(-1,17]`,"before 18-40"=`(17,40]`,"before >40"=`(40,120]`) %>% 
  left_join(tmp2 %>% 
              filter(period == "before") %>% 
              select(dx,dx_ver,"before All Ages"=frac)) %>% 
  arrange(desc(`before All Ages`))

during_counts %>% 
  left_join(before_counts) %>% 
  arrange(desc(`All Ages`)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/results/aug_2025/ssd_counts_by_age.csv")


##################
#### DX Plots ####
##################

plot_path <- "~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/results/aug_2025/"

### ALL DX ---------------------------------------------------------------------

tmp <- dx_week_counts %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

pdf(paste0(plot_path,"all_ssds_by_age.pdf"),onefile = TRUE)
for (i in 1:ceiling(nrow(tmp)/8)){
  
  tmp_conds <- tmp %>%
    slice((i*8-7):(i*8)) %>% 
    select(dx,dx_ver)
  
  p <- dx_week_counts %>%
    inner_join(tmp_conds, by = c("dx","dx_ver")) %>% 
    ggplot(aes(week,frac,color = age_cat)) +
    geom_line() +
    facet_wrap(~label, nrow = 4,scales = "free_y") +
    scale_x_reverse() +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("% of patients with visit") +
    xlab("weeks before diagnosis")
  
  print(p)
}
dev.off()


### DX-9 -----------------------------------------------------------------------

tmp <- dx_week_counts %>% 
  filter(dx_ver==9) %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

pdf(paste0(plot_path,"icd9_ssds_by_age.pdf"),onefile = TRUE)
for (i in 1:ceiling(nrow(tmp)/8)){
  
  tmp_conds <- tmp %>%
    slice((i*8-7):(i*8)) %>% 
    select(dx,dx_ver)
  
  p <- dx_week_counts %>%
    inner_join(tmp_conds, by = c("dx","dx_ver")) %>% 
    ggplot(aes(week,frac,color = age_cat)) +
    geom_line() +
    facet_wrap(~label, nrow = 4,scales = "free_y") +
    scale_x_reverse() +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("% of patients with visit") +
    xlab("weeks before diagnosis")
  
  print(p)
}
dev.off()


### DX-10 ----------------------------------------------------------------------

tmp <- dx_week_counts %>% 
  filter(dx_ver==10) %>% 
  group_by(dx,dx_ver) %>% 
  summarise(n = sum(n)) %>% 
  arrange(desc(n)) %>% 
  ungroup()

pdf(paste0(plot_path,"icd10_ssds_by_age.pdf"),onefile = TRUE)
for (i in 1:ceiling(nrow(tmp)/8)){
  
  tmp_conds <- tmp %>%
    slice((i*8-7):(i*8)) %>% 
    select(dx,dx_ver)
  
  p <- dx_week_counts %>%
    inner_join(tmp_conds, by = c("dx","dx_ver")) %>% 
    ggplot(aes(week,frac,color = age_cat)) +
    geom_line() +
    facet_wrap(~label, nrow = 4,scales = "free_y") +
    scale_x_reverse() +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("% of patients with visit") +
    xlab("weeks before diagnosis")
  
  print(p)
}
dev.off()


