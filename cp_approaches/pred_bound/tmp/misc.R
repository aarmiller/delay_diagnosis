library(tidyverse)

rm(list = ls())

load("/Shared/AML/params/final_delay_params.RData")

cond_name <- "sarcoid_lung"

delay_params <- final_delay_params[[cond_name]]


## Load Data -------------------------------------------------------------------
# load index cases
load(paste0(delay_params$out_path,"index_cases.RData"))

# extract patient ids and number of patients
patient_ids <- index_cases %>% 
  distinct(patient_id)

data_in_path <- paste0(delay_params$base_path,"delay_results/")
load(paste0(data_in_path,"all_dx_visits.RData"))

load(paste0(data_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes("sarcoid_lung") %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver) %>% 
  distinct() %>% 
  filter(!is.na(dx))

all_dx_visits <- inner_join(all_dx_visits, patient_ids, by = join_by("patient_id"))
ssd_dx_visits <- inner_join(all_dx_visits, ssd_codes, by = c("dx", "dx_ver"))


#### Setup outer bootstrap data ------------------------------------------------

# sample overall set of patients for each bootstrap
set.seed(192837)
tmp <- tibble(boot_trial=1:delay_params$boot_trials) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(patient_id = sample(patient_ids$patient_id, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# tmp$boot_sample[[1]]


#### compute counts for ssd visits ---------------------------------------------
get_ssd_counts <- function(patient_set){
  patient_set %>% 
    inner_join(ssd_dx_visits, by = join_by("patient_id"),relationship = "many-to-many") %>%
    mutate(period = -days_since_index) %>%
    distinct(patient_id,boot_id,period) %>% 
    count(period) %>% 
    filter(period>0) %>% 
    mutate(dow = as.factor(period %% 7))
}

get_ssd_counts(tmp$boot_sample[[1]])

# for each set of patients pull the ssd counts
tmp <- tmp %>% 
  mutate(ssd_vis_count = map(boot_sample,get_ssd_counts))

# tmp$all_vis_count[[1]]
# tmp$ssd_vis_count[[1]]

save(tmp,file = "/Shared/AML/tmp_transfer/saroid_lung_bootstrap.RData")

# rm(list = ls())
load("/Volumes/AML/tmp_transfer/saroid_lung_bootstrap.RData")

tmp %>% 
  select(boot_trial,ssd_vis_count) %>% 
  unnest(ssd_vis_count) %>% 
  mutate(boot_trial==as.factor(boot_trial)) %>% 
  # filter(boot_trial==5) %>% 
  ggplot(aes(period,n,group = boot_trial)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  scale_x_reverse()


tmp$ssd_vis_count[[1]]

### fit models -----------------------------------------------------------------
fit_models <- function(count_data){
  models %>% 
    mutate(counts = map2(model,cp,
                         ~return_fits(data = count_data,
                                      model = .x,
                                      cp = .y,
                                      periodicity = delay_params$periodicity)))
}

# fit trends for both all visits and ssd visits
tmp <- tmp %>% 
  mutate(all_vis_count = map(all_vis_count,fit_models)) %>% 
  mutate(ssd_vis_count = map(ssd_vis_count,fit_models))

boot_data <- tmp