
rm(list=ls()[!(ls() %in% c("proj_name", "project_name"))])
library(tidyverse)

# args = commandArgs(trailingOnly=TRUE)

# name of condition
# proj_name <- args[1]
# proj_name <-  "dengue"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

### Load delay params ----------------------------------------------------------

load("/Shared/AML/params/final_delay_params_kaiser.RData")

delay_params <- final_delay_params[[proj_name]]

out_path <- delay_params$out_path

if (!dir.exists(out_path)){
  dir.create(out_path)
}else{
  print("dir exists")
}

### Paths for plot output ------------------------------------------------------

plot_path <- paste0(out_path,"potential_ssd_plots/")

if (!dir.exists(plot_path)){
  dir.create(plot_path)
}else{
  print("dir exists")
}

### Load Data ------------------------------------------------------------------
load(paste0(delay_params$out_path,"index_cases.RData"))

db_path <- delay_params$small_db_path

db <- src_sqlite(list.files(db_path, full.names = T))

patient_ids <- index_cases %>% distinct(patient_id)

index_dx_dates <- db %>% tbl("index_dx_dates") %>% 
  filter(patient_id %in% local(patient_ids$patient_id)) %>% collect() %>% 
  select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id") 
  
all_dx_visits <- db %>% tbl("all_dx_visits") %>% 
  filter(patient_id %in% local(patient_ids$patient_id)) %>% collect() %>% 
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) 

# subset to visits that did not come from other encounter type
load(paste0(delay_base_path,"delay_tm.RData"))
tm <- tm %>%   
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>%
  mutate(days_since_index = days_since_index - shift) %>%
  select(-shift) %>%
  # filter(days_since_index<=0) %>% # we want visit beyond index here
  select(patient_id, days_since_index, outpatient, ed, inpatient, other) %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0)) # subset to only visits from AV, ED, IP, or IS

all_dx_visits <- all_dx_visits %>% 
  inner_join(tm %>% distinct(patient_id, days_since_index), by = c("patient_id", "days_since_index")) # visit days from other encounter types only removed

all_dx_visits <- all_dx_visits %>% inner_join(patient_ids)
index_dx_dates <- index_dx_dates%>% inner_join(patient_ids)


index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=delay_params$upper_bound)

tmp_dx_visits <- all_dx_visits %>%
  filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
  inner_join(distinct(index_dx_dates, patient_id))

tmp_dx_visits9 <- tmp_dx_visits %>% filter(dx_ver=="09")
tmp_dx_visits10 <- tmp_dx_visits %>% filter(dx_ver=="10")
tmp_dx_visitsOT <- tmp_dx_visits %>% filter(dx_ver=="OT")

### Make All visits plot -------------------------------------------------------
tmp_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  geom_smooth( color = "red",span = 0.25) +
  theme_minimal() +
  ggtitle(paste0("All visits before ",cond_name))
ggsave(paste0(plot_path,"all_visits_before.pdf"),width = 6, height = 5)


### Get counts of top codes ----------------------------------------------------

dx9_counts <- tmp_dx_visits9 %>%
  distinct(patient_id,dx,days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  count(dx) %>%
  arrange(desc(n)) %>%
  left_join(tibble(icd::icd9cm_hierarchy) %>%
              select(dx=code,desc=long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  mutate(desc = paste0(dx," - ",str_sub(desc,start = 0,end = 40)))

dx9_counts %>% 
  slice(1:1000) %>% 
  write_csv(paste0(plot_path,"icd9_top1000.csv"))


dx10_counts <- tmp_dx_visits10 %>%
  distinct(patient_id,dx,days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  count(dx) %>%
  arrange(desc(n)) %>%
  left_join(tibble(icd::icd10cm2019) %>%
              select(dx=code,desc=long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  mutate(desc = paste0(dx," - ",str_sub(desc,start = 0,end = 40)))

dx10_counts %>% 
  slice(1:1000) %>% 
  write_csv(paste0(plot_path,"icd10_top1000.csv"))

dxOT_counts <- tmp_dx_visitsOT %>%
  distinct(patient_id,dx,days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  count(dx) %>%
  arrange(desc(n)) %>%
  left_join(tibble(icd::icd10cm2019) %>%
              select(dx=code,desc=long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  mutate(desc = paste0(dx," - ",str_sub(desc,start = 0,end = 40)))

dxOT_counts %>% 
  slice(1:1000) %>% 
  write_csv(paste0(plot_path,"icdOT_top1000.csv"))

### Plot dx 9  before ----------------------------------------------------------

pdf(paste0(plot_path,"dx9_before_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dx9_counts)/8))){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))
  
  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) 
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n)) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      theme_minimal()
  }
  
  print(p)
}
dev.off()


### Plot dx 10  before ----------------------------------------------------------

pdf(paste0(plot_path,"dx10_before_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dx10_counts)/8))){
  tmp_conds <- dx10_counts %>%
    slice((i*8-7):(i*8))
  
  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits10, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) 
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n)) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      theme_minimal()
  }
  
  print(p)
}
dev.off()


### Plot dx OT  before ----------------------------------------------------------

pdf(paste0(plot_path,"dxOT_before_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dxOT_counts)/8))){
  tmp_conds <- dxOT_counts %>% filter(n >1) %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visitsOT, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) 
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n)) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      theme_minimal()
  }
  
  print(p)
}
dev.off()


### Plot dx 9  before and after ------------------------------------------------

pdf(paste0(plot_path,"dx9_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dx9_counts)/8))){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))
  
  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After"))
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n, color = period)) +
      geom_line() +
      geom_smooth(aes(group = period),color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      geom_vline(aes(xintercept = 0)) +
      theme_minimal() +
      theme(legend.position="none")
  }
  
  print(p)
}
dev.off()


### Plot dx 10  before and after ------------------------------------------------

pdf(paste0(plot_path,"dx10_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dx10_counts)/8))){
  tmp_conds <- dx10_counts %>%
    slice((i*8-7):(i*8))
  
  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits10, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) 
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n, color = period)) +
      geom_line() +
      geom_smooth(aes(group = period),color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      geom_vline(aes(xintercept = 0)) +
      theme_minimal() +
      theme(legend.position="none")
  }

  print(p)
}
dev.off()

### Plot dx OT  before and after ------------------------------------------------

pdf(paste0(plot_path,"dxOT_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:(ceiling(nrow(dxOT_counts)/8))){
  tmp_conds <- dxOT_counts %>%
    slice((i*8-7):(i*8))
  
  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visitsOT, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) 
  
  if(length(unique(p$desc))>1){
    p <- p %>%
      ggplot(aes(days_since_index,n, color = period)) +
      geom_line() +
      geom_smooth(aes(group = period),color = "red",span = 0.25) +
      facet_wrap(~desc, scales = "free_y", nrow =4) +
      geom_vline(aes(xintercept = 0)) +
      theme_minimal() +
      theme(legend.position="none")
  }
  
  print(p)
}
dev.off()
