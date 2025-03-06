
rm(list = ls())

library(tidyverse)


pd_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pd/pd.db")
als_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/als/als.db")
mg_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/mg/mg.db")

ssd_codes <- c("7820", "7804", "28521", "5789", "2767", "29620", 
                           "5849", "E8889", "30000", "79902", "78451", "78096",
                           "7802", "7282", "78079", "27651", "72989", "2809",
                           "49121", "78459", "V573", "2761", "V1588", "58881",
                           "5859", "25062", "7813", "40391", "5856", "78052",
                           "2760", "4010", "4019", "40390", "7843", "78820", 
                           "5990", "78060", "78720", "2768")

pd_ssd_vis <- pd_db %>% 
  tbl("all_dx_visits") %>% 
  filter(dx %in% ssd_codes) %>% 
  collect()

als_ssd_vis <- als_db %>% 
  tbl("all_dx_visits") %>% 
  filter(dx %in% ssd_codes) %>% 
  collect()

mg_ssd_vis <- mg_db %>% 
  tbl("all_dx_visits") %>% 
  filter(dx %in% ssd_codes) %>% 
  collect()

pd_index <- pd_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

als_index <- als_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

mg_index <- mg_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

save(pd_ssd_vis,als_ssd_vis,mg_ssd_vis,
     pd_index,als_index,mg_index,
     file = "/Shared/AML/tmp_transfer/jacob_pd_grant.RData")


### Run Locally --------

load("/Volumes/AML/tmp_transfer/jacob_pd_grant.RData")


all_ssd_vis <- bind_rows(mutate(pd_ssd_vis,group = "PD"),
          mutate(als_ssd_vis,group = "ALS"),
          mutate(mg_ssd_vis,group = "MG"))

ssd_counts <- all_ssd_vis %>% 
  distinct(patient_id,group,dx,days_since_index) %>% 
  count(group,days_since_index,dx) %>% 
  filter(between(days_since_index,-365,-1))


any_ssd_counts <- all_ssd_vis %>% 
  distinct(patient_id,group,days_since_index) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  count(group,days_since_index)


any_ssd_counts %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_wrap(~group, scales = "free_y")


all_index <- bind_rows(mutate(pd_index,group = "PD"),
                         mutate(als_index,group = "ALS"),
                         mutate(mg_index,group = "MG"))

pop_counts <- all_index %>% 
  distinct(group, patient_id) %>% 
  count(group,name = "pop_n")


save(all_ssd_vis,all_index)

## Plot SSD incidence ------
any_ssd_counts %>% 
  inner_join(pop_counts) %>% 
  mutate(incidence = 1000*n/pop_n) %>% 
  ggplot(aes(days_since_index,incidence, color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth() +
  theme_minimal()


######################
#### Any SSD Code ####
######################

### Filter to different enrollment period ---------

obs_window <- 2*365

plot_any_ssd <- function(obs_window){
  include_ids <- all_index %>% 
    filter(time_before_index>=obs_window) %>% 
    distinct(patient_id,group)
  
  include_pop_count <- include_ids %>% 
    count(group, name = "pop_n")
  
  all_ssd_vis %>% 
    inner_join(include_ids,by = join_by(patient_id, group)) %>% 
    filter(between(days_since_index,-obs_window,-1)) %>% 
    distinct(patient_id,days_since_index,group) %>% 
    count(days_since_index,group) %>% 
    inner_join(include_pop_count,by = join_by(group)) %>% 
    mutate(incidence = 1000*n/pop_n) %>% 
    ggplot(aes(-days_since_index,incidence, color = group)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth() +
    theme_minimal() +
    scale_x_reverse() +
    ylab("Daily Visit Incidence (per 1000K enrollees)") +
    xlab("Days Before Diagnosis") +
    theme(legend.title = element_blank()) +
    ggtitle(paste0("SSD Visits Incidence: 1 to ",obs_window," days before diagnosis"))
    
}

plot_any_ssd(365)
plot_any_ssd(2*365)
plot_any_ssd(3*365)


make_any_ssd_window_counts <- function(obs_window){
  
  include_ids <- all_index %>% 
    filter(time_before_index>=obs_window) %>% 
    distinct(patient_id,group)
  
  include_pop_count <- include_ids %>% 
    count(group, name = "pop_n")
  
  all_ssd_vis %>% 
    inner_join(include_ids,by = join_by(patient_id, group)) %>% 
    filter(between(days_since_index,-obs_window,-1)) %>% 
    distinct(patient_id,days_since_index,group) %>% 
    count(group,days_since_index) %>% 
    inner_join(include_pop_count, by = "group") %>% 
    mutate(incidence = 1000*n/pop_n)
  
}

any_ssd_counts <- tibble(obs_window = 365*c(1,2,3,4)) %>% 
  mutate(ssd_counts = map(obs_window,make_ssd_window_counts))

any_ssd_counts <- any_ssd_counts %>% 
  unnest(ssd_counts)

any_ssd_counts %>% 
  mutate(obs_window = obs_window/365) %>% 
  mutate(obs_window = paste0(obs_window," Years")) %>% 
  mutate(obs_window = fct_relevel(obs_window, c("1 Years",
                                                "2 Years",
                                                "3 Years",
                                                "4 Years"))) %>% 
  ggplot(aes(-days_since_index,incidence, color = group)) +
  geom_point(alpha = 0.25, size = 0.5) +
  geom_smooth(size = 1.5) +
  theme_minimal() +
  scale_x_reverse() +
  ylab("Daily Visit Incidence (per 1000K enrollees)") +
  xlab("Days Before Diagnosis") +
  theme(legend.title = element_blank()) +
  facet_wrap(~obs_window, scale = "free")


#### By SSD ####


make_ssd_window_counts <- function(obs_window){
  
  include_ids <- all_index %>% 
    filter(time_before_index>=obs_window) %>% 
    distinct(patient_id,group)
  
  include_pop_count <- include_ids %>% 
    count(group, name = "pop_n")
  
  all_ssd_vis %>% 
    inner_join(include_ids,by = join_by(patient_id, group)) %>% 
    filter(between(days_since_index,-obs_window,-1)) %>% 
    distinct(patient_id,dx,days_since_index,group) %>% 
    count(group,dx,days_since_index) %>% 
    inner_join(include_pop_count, by = "group") %>% 
    mutate(incidence = 1000*n/pop_n)
  
}

rm(obs_window,include_ids,include_pop_count)

individual_ssd_counts <- tibble(obs_window = 365*c(1,1.5,2,2.5,3,3.5,4)) %>% 
  mutate(ssd_counts = map(obs_window,make_ssd_window_counts))

individual_ssd_counts %>% 
  unnest(ssd_counts) %>% 
  filter(obs_window == 3*365) %>% 
  ggplot(aes(-days_since_index,incidence, color = group)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth() +
  theme_minimal() +
  scale_x_reverse() +
  ylab("Daily Visit Incidence (per 1000K enrollees)") +
  xlab("Days Before Diagnosis") +
  theme(legend.title = element_blank()) +
  facet_wrap(~dx, scales = "free_y") +
  



individual_ssd_counts %>% 
  unnest(ssd_counts) %>% 
  filter(obs_window == 3*365) %>% 
  filter(dx %in% c("2809","28521","4019","5856","58881,5990",
                   "7802","7804")) %>% 
  inner_join(filter(codeBuildr::all_icd_labels,dx_ver==9)) %>% 
  mutate(desc = paste0(dx,": ", desc)) %>% 
  ggplot(aes(-days_since_index,incidence, color = group)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth() +
  theme_minimal() +
  scale_x_reverse() +
  ylab("Daily Visit Incidence (per 1000K enrollees)") +
  xlab("Days Before Diagnosis") +
  theme(legend.title = element_blank()) +
  facet_wrap(~desc, scales = "free_y") +
  theme(strip.text = element_text(size = 12))




