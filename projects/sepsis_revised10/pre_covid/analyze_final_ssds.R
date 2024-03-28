

library(tidyverse)
library(icd)
library(codeBuildr)


cond_name <- "sepsis_pre_covid"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)

sim_res_path <- paste0(delay_params$out_path,"sim_results/exponential_cp14/")


###################
#### Load Data ####
###################


## load in dx visits -----------------------------------------------------------
load(paste0(delay_params$base_path,"delay_results/all_dx_visits.RData"))




## load in timemap -------------------------------------------------------------
# load(paste0(delay_params$base_path,"delay_results/delay_tm.RData"))
# # tm
# tm <- tm %>% inner_join(select(index_cases,patient_id), by = "patient_id")

## load simulation results -----------------------------------------------------
load(paste0(sim_res_path,"sim_res_ssd.RData"))
# sim_res_ssd


########################
#### All SSD Counts ####
########################

ssd_codes <- codeBuildr::load_ssd_codes("sepsis_revised10") %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver) %>%
  distinct() %>%
  filter(!is.na(dx))


count_ssds <- function(data){
  data %>% 
    inner_join(sim_obs,by = join_by(obs)) %>%
    inner_join(ssd_dx_visits,by = join_by(patient_id, days_since_index),relationship = "many-to-many") %>% 
    distinct(dx,boot_id,days_since_index) %>% 
    count(dx) %>% 
    arrange(desc(n)) %>% 
    left_join(filter(codeBuildr::all_icd_labels,dx_ver==10),by = join_by(dx))
}

count_ssds(sim_res_ssd$res[[1]])

tmp <- sim_res_ssd %>% 
  mutate(ssd_count = map(res,count_ssds))

ssd_counts <- tmp

ssd_counts <- ssd_counts %>% 
  mutate(miss_count = map_int(res,nrow)) %>% 
  select(sim_trial,boot_trial,ssd_count,miss_count) %>% 
  unnest(ssd_count)

ssd_counts_summary <- ssd_counts %>% 
  mutate(miss_frac = 100*n/miss_count) %>% 
  group_by(dx,desc) %>% 
  summarise(mean_miss_frac = mean(miss_frac),
            low_miss_frac = quantile(miss_frac,probs =0.025),
            high_miss_frac = quantile(miss_frac,probs =0.975),
            median_miss_frac = median(miss_frac),
            min_miss_frac = min(miss_frac),
            max_miss_frac = max(miss_frac)) %>% 
  arrange(desc(mean_miss_frac)) %>% 
  ungroup()

save(ssd_counts_summary,file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts.RData")

load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts.RData")

ssd_counts_summary %>% 
  mutate_at(vars(mean_miss_frac:max_miss_frac),~round(.,2)) %>% 
  mutate(out = paste0(mean_miss_frac," (",low_miss_frac,"-",high_miss_frac,")")) %>% 
  select(`ICD Code`=dx,Description = desc, `Frequency (95% CI)`=out) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/results/ssd_freq_count.csv")

ssd_counts_summary %>% 
  summarise(sum(mean_miss_frac))


##############################
##### Using Custom Labels ####
##############################

### Category 1 -----------------------------------------------------------------

# get aggregated SSD categories
ssd_cats <- read_csv("/Shared/AML/cdc_sepsis/delay_dx/data/sepsis_ssd_categories.csv")

ssd_cats_corretion <- read_csv("/Shared/AML/cdc_sepsis/delay_dx/data/corrected_labels.csv")

ssd_codes <- ssd_cats_corretion %>%
  rename(group1=label) %>% 
  inner_join(ssd_cats) %>% 
  select(dx=`ICD Code`,label=`new lab`) %>% 
  mutate(dx_ver=10)


ssd_dx_visits <- inner_join(all_dx_visits, ssd_codes, by = c("dx", "dx_ver"))


count_ssd_cats <- function(data){
  data %>% 
    inner_join(sim_obs,by = join_by(obs)) %>%
    inner_join(ssd_dx_visits,by = join_by(patient_id, days_since_index),relationship = "many-to-many") %>% 
    distinct(label,boot_id,days_since_index) %>% 
    count(label) %>% 
    arrange(desc(n)) %>% 
    mutate(frac = n/nrow(data))
}


count_ssd_cats(sim_res_ssd$res[[1]])

tmp <- sim_res_ssd %>% 
  mutate(ssd_count = map(res,count_ssd_cats))


ssd_counts <- tmp %>% 
  select(sim_trial,boot_trial,ssd_count) %>% 
  unnest(ssd_count)

ssd_counts_summary <- ssd_counts %>% 
  mutate(miss_frac = 100*frac) %>% 
  group_by(label) %>% 
  summarise(mean_miss_frac = mean(miss_frac),
            low_miss_frac = quantile(miss_frac,probs =0.025),
            high_miss_frac = quantile(miss_frac,probs =0.975),
            median_miss_frac = median(miss_frac),
            min_miss_frac = min(miss_frac),
            max_miss_frac = max(miss_frac)) %>% 
  arrange(desc(mean_miss_frac)) %>% 
  ungroup()

save(ssd_counts_summary,file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_cat1.RData")


### Category 2 -----------------------------------------------------------------


# get aggregated SSD categories
ssd_cats <- read_csv("/Shared/AML/cdc_sepsis/delay_dx/data/sepsis_ssd_categories.csv")

ssd_codes <- ssd_cats %>% 
  select(dx=`ICD Code`,label=group2) %>% 
  mutate(dx_ver=10)


ssd_dx_visits <- inner_join(all_dx_visits, ssd_codes, by = c("dx", "dx_ver"))


count_ssd_cats <- function(data){
  data %>% 
    inner_join(sim_obs,by = join_by(obs)) %>%
    inner_join(ssd_dx_visits,by = join_by(patient_id, days_since_index),relationship = "many-to-many") %>% 
    distinct(label,boot_id,days_since_index) %>% 
    count(label) %>% 
    arrange(desc(n)) %>% 
    mutate(frac = n/nrow(data))
}


count_ssd_cats(sim_res_ssd$res[[1]])

tmp <- sim_res_ssd %>% 
  mutate(ssd_count = map(res,count_ssd_cats))


ssd_counts <- tmp %>% 
  select(sim_trial,boot_trial,ssd_count) %>% 
  unnest(ssd_count)

ssd_counts_summary <- ssd_counts %>% 
  mutate(miss_frac = 100*frac) %>% 
  group_by(label) %>% 
  summarise(mean_miss_frac = mean(miss_frac),
            low_miss_frac = quantile(miss_frac,probs =0.025),
            high_miss_frac = quantile(miss_frac,probs =0.975),
            median_miss_frac = median(miss_frac),
            min_miss_frac = min(miss_frac),
            max_miss_frac = max(miss_frac)) %>% 
  arrange(desc(mean_miss_frac)) %>% 
  ungroup()

save(ssd_counts_summary,file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_cat2.RData")


### Pain -----------------------------------------------------------------------


# get aggregated SSD categories
ssd_cats <- read_csv("/Shared/AML/cdc_sepsis/delay_dx/data/sepsis_ssd_categories.csv")

ssd_codes <- ssd_cats %>% 
  filter(`Pain Yes/No`=="yes") %>% 
  select(dx=`ICD Code`) %>% 
  mutate(label = "pain") %>% 
  mutate(dx_ver=10)


ssd_dx_visits <- inner_join(all_dx_visits, ssd_codes, by = c("dx", "dx_ver"))


count_ssd_cats <- function(data){
  data %>% 
    inner_join(sim_obs,by = join_by(obs)) %>%
    inner_join(ssd_dx_visits,by = join_by(patient_id, days_since_index),relationship = "many-to-many") %>% 
    distinct(label,boot_id,days_since_index) %>% 
    count(label) %>% 
    arrange(desc(n)) %>% 
    mutate(frac = n/nrow(data))
}

tmp <- sim_res_ssd %>% 
  mutate(ssd_count = map(res,count_ssd_cats))


ssd_counts <- tmp %>% 
  select(sim_trial,boot_trial,ssd_count) %>% 
  unnest(ssd_count)

ssd_counts_summary <- ssd_counts %>% 
  mutate(miss_frac = 100*frac) %>% 
  group_by(label) %>% 
  summarise(mean_miss_frac = mean(miss_frac),
            low_miss_frac = quantile(miss_frac,probs =0.025),
            high_miss_frac = quantile(miss_frac,probs =0.975),
            median_miss_frac = median(miss_frac),
            min_miss_frac = min(miss_frac),
            max_miss_frac = max(miss_frac)) %>% 
  arrange(desc(mean_miss_frac)) %>% 
  ungroup()

save(ssd_counts_summary,file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_pain.RData")

rm(list = ls())

#######################
#### Export Tables ####
#######################

load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_cat1.RData")
cat1_res <- ssd_counts_summary
load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_cat2.RData")
cat2_res <- ssd_counts_summary
load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/ssd_summary_counts_pain.RData")
pain_res <- ssd_counts_summary

tab1 <- bind_rows(cat1_res,pain_res) %>% 
  arrange(desc(median_miss_frac)) %>% 
  mutate(res =paste0(round(median_miss_frac,2)," (",round(low_miss_frac,2),"-",round(high_miss_frac,2),")")) %>% 
  select(label,res)

tab1 %>% write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/manuscripts/ssd_res_table1.csv")

tab1 %>% arrange(label) %>% write_csv("~/Downloads/tmp.csv")


tab2 <- bind_rows(cat2_res,pain_res) %>% 
  arrange(desc(median_miss_frac)) %>% 
  mutate(res =paste0(round(median_miss_frac,2)," (",round(low_miss_frac,2),"-",round(high_miss_frac,2),")")) %>% 
  select(label,res)

tab2 %>% write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/manuscripts/ssd_res_table2.csv")


## Appendix table with codes

ssd_cats <- read_csv("/Volumes/AML/cdc_sepsis/delay_dx/data/sepsis_ssd_categories.csv")

ssd_cats_corretion <- read_csv("/Volumes/AML/cdc_sepsis/delay_dx/data/corrected_labels.csv")

ssd_cats_corretion
tmp1 <- ssd_cats %>% 
  select(label = group1,dx = `ICD Code`) %>% 
  inner_join(ssd_cats_corretion) %>% 
  select(label = `new lab`,dx) %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(codes = map_chr(data,~.$dx %>% paste0(collapse=", "))) %>% 
  select(label,codes)

tmp2 <- tibble(label = "pain",
       codes = ssd_cats %>% 
         filter(`Pain Yes/No`=="yes") %>% 
         .$`ICD Code` %>% 
         paste0(collapse = ", "))

bind_rows(tmp1,tmp2) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/manuscripts/1_19_24/code_table.csv")
