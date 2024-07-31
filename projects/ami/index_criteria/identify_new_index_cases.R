

library(tidyverse)


## Pull in original timemap ----------------------------------------------------

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/ami/ami.db")

tm <- db %>% tbl("tm") %>% collect()

## Pull in all visit information -----------------------------------------------
load("/Shared/Statepi_Diagnosis/prelim_results/ami/delay_results/all_dx_visits.RData")

# merge in index and filter to visits within 14 days of diagnosis
tm <- select(index_dx_dates,patient_id,index_date) %>% 
  inner_join(tm) %>% 
  mutate(days_since_index = svcdate-index_date) %>% 
  filter(between(days_since_index,-14,14))

## Identify two cases meeting inclusion criteria -------------------------------

# cases diagnosed on inpatient date
index_pat1 <- tm %>% 
  filter(days_since_index==0 & inpatient==1) %>% 
  distinct(patient_id)

# cases diagnosed within 7 days of inpatient admission
index_pat2 <- tm %>% 
  filter(days_since_index==0 & inpatient==0) %>% 
  distinct(patient_id) %>% 
  inner_join(tm) %>% 
  filter(between(days_since_index,0,7) & inpatient==1) %>% 
  distinct(patient_id)

# combine and add indicator for how identified (may want to look at this later)
final_patients <- bind_rows(mutate(index_pat1,group = 1L),mutate(index_pat2,group = 2)) 


## Split patients into male and female -----------------------------------------
load("/Shared/Statepi_Diagnosis/prelim_results/ami/delay_results/demo_data.RData")

final_patients <- final_patients %>% 
  inner_join(select(demo1,patient_id,sex))


save(final_patients, file = "/Shared/Statepi_Diagnosis/projects/ami/final_patients.RData")



## Compute SSD counts for each group -------------------------------------------

ssd_codes <- codeBuildr::load_ssd_codes("ami") %>% 
  mutate(dx = code, dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
  select(dx,dx_ver)

tmp <- all_dx_visits %>% 
  inner_join(ssd_codes) %>% 
  distinct(patient_id,days_since_index)

sex_ssd_counts <- tmp %>% 
  inner_join(final_patients) %>% 
  count(sex,days_since_index)

save(sex_ssd_counts, file = "/Shared/Statepi_Diagnosis/projects/ami/tmp/sex_ssd_counts.RData")


load("/Volumes/Statepi_Diagnosis/projects/ami/tmp/sex_ssd_counts.RData")

save(sex_ssd_counts, file = "~/Data/Statepi_Diagnosis/projects/ami/sex_ssd_counts.RData")


sex_ssd_counts %>% 
  mutate(sex = ifelse(sex == 1, "Male","Female")) %>% 
  filter(days_since_index< -1) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  scale_x_reverse() +
  facet_wrap(~sex, scale = "free_y") +
  theme_bw() +
  geom_vline(aes(xintercept = 60), linetype = 2) +
  xlab("Days Before Index")


