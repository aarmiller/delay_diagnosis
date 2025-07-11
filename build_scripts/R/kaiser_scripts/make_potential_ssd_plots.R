rm(list=ls()[!(ls() %in% c("cond_name"))])
library(tidyverse)

### Load delay params ----------------------------------------------------------

# Load Delay Params
load("/Shared/AML/params/delay_any_params_kaiser.RData")

delay_params <- delay_any_params[[cond_name]]

base_path <- delay_params$path

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name)

if (!dir.exists(out_path)){
  dir.create(out_path)
}else{
  print("dir exists")
}


### Paths for plot output ------------------------------------------------------

plot_path <- paste0(out_path,"/potential_ssd_plots/")

if (!dir.exists(plot_path)){
  dir.create(plot_path)
}else{
  print("dir exists")
}

### Load Data ------------------------------------------------------------------

db <- src_sqlite(paste0(base_path,"/",cond_name,".db"))

index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()
all_dx_visits <- db %>% tbl("all_dx_visits") %>% collect()

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
for (i in 1:125){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    geom_smooth( color = "red",span = 0.25) +
    facet_wrap(~desc, scales = "free_y", nrow =4) +
    theme_minimal()

  print(p)
}
dev.off()


### Plot dx 10  before ----------------------------------------------------------

pdf(paste0(plot_path,"dx10_before_top_1000.pdf"),onefile = TRUE)
for (i in 1:125){
  tmp_conds <- dx10_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits10, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    geom_smooth( color = "red",span = 0.25) +
    facet_wrap(~desc, scales = "free_y", nrow =4) +
    theme_minimal()

  print(p)
}
dev.off()


### Plot dx OT  before ----------------------------------------------------------

# pdf(paste0(plot_path,"dxOT_before_top_1000.pdf"),onefile = TRUE)
# for (i in 1:125){
#   tmp_conds <- dxOT_counts %>% filter(n >1) %>% 
#     slice((i*8-7):(i*8))
#   
#   p <- tmp_conds %>%
#     select(dx,desc) %>%
#     inner_join(tmp_dx_visitsOT, by ="dx") %>%
#     distinct(patient_id,dx,desc,days_since_index) %>%
#     count(dx,desc,days_since_index) %>%
#     filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
#     ggplot(aes(days_since_index,n)) +
#     geom_line() +
#     geom_smooth( color = "red",span = 0.25) +
#     facet_wrap(~desc, scales = "free_y", nrow =4) +
#     theme_minimal()
#   
#   print(p)
# }
# dev.off()


### Plot dx 9  before and after ------------------------------------------------

pdf(paste0(plot_path,"dx9_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:125){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,n, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~desc, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none")

  print(p)
}
dev.off()


### Plot dx 10  before and after ------------------------------------------------

pdf(paste0(plot_path,"dx10_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:125){
  tmp_conds <- dx10_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits10, by ="dx") %>%
    distinct(patient_id,dx,desc,days_since_index) %>%
    count(dx,desc,days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,n, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~desc, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none")

  print(p)
}
dev.off()

