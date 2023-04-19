
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "lung_cancer"

### Load delay params ----------------------------------------------------------

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

base_path <- delay_params$path


### Paths for plot output ------------------------------------------------------

plot_path <- paste0(base_path,"/potential_ssd_plots/")

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
  inner_join(distinct(index_dx_dates, enrolid))

tmp_dx_visits9 <- tmp_dx_visits %>% filter(dx_ver==9)
tmp_dx_visits10 <- tmp_dx_visits %>% filter(dx_ver==10)

### Make All visits plot -------------------------------------------------------
tmp_dx_visits %>%
  distinct(enrolid,days_since_index) %>%
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
  distinct(enrolid,dx,days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  count(dx) %>%
  arrange(desc(n)) %>%
  left_join(tibble(icd::icd9cm_hierarchy) %>%
              select(dx=code,desc=long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  mutate(desc = paste0(dx," - ",str_sub(desc,start = 0,end = 40)))


dx10_counts <- tmp_dx_visits10 %>%
  distinct(enrolid,dx,days_since_index) %>%
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
  count(dx) %>%
  arrange(desc(n)) %>%
  left_join(tibble(icd::icd10cm2019) %>%
              select(dx=code,desc=long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  mutate(desc = paste0(dx," - ",str_sub(desc,start = 0,end = 40)))

### Plot dx 9  before ----------------------------------------------------------

pdf(paste0(plot_path,"dx9_before_top_1000.pdf"),onefile = TRUE)
for (i in 1:125){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(enrolid,dx,desc,days_since_index) %>%
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
    distinct(enrolid,dx,desc,days_since_index) %>%
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


### Plot dx 9  before and after ------------------------------------------------

pdf(paste0(plot_path,"dx9_before_after_top_1000.pdf"),onefile = TRUE)
for (i in 1:125){
  tmp_conds <- dx9_counts %>%
    slice((i*8-7):(i*8))

  p <- tmp_conds %>%
    select(dx,desc) %>%
    inner_join(tmp_dx_visits9, by ="dx") %>%
    distinct(enrolid,dx,desc,days_since_index) %>%
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
    distinct(enrolid,dx,desc,days_since_index) %>%
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

