
rm(list=ls()[!(ls() %in% c("cond_name"))])
library(tidyverse)
library(bit64)
# library(parallel)
library(smallDB)
# library(delaySim)
library(lubridate)

############################################
#### Load Data and Prepare Output paths ####
############################################

# Load Delay Params
load("/Shared/AML/params/delay_any_params_kaiser.RData")

delay_params <- delay_any_params[[cond_name]]

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name)
# out_path <- paste0("~/Data/tmp/",cond_name)

if (!dir.exists(out_path)){
  dir.create(out_path)
}

# connect to database
con <- DBI::dbConnect(RSQLite::SQLite(), paste0(delay_params$path,cond_name,".db"))
# con <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

#### Add output folders ---------------------------------------------------------

sim_out_path <- paste0(out_path,"/delay_results")

if (!dir.exists(sim_out_path)){
  dir.create(sim_out_path)
}

cluster_out_path <- paste0(out_path,"/cluster_results")

if (!dir.exists(cluster_out_path)){
  dir.create(cluster_out_path)
}

##################################
##### Get Visit Data Extracts ####
##################################

#### Index DX dates ------------------------------------------------------------

# Filter to patient_ids that can be observed 
index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=delay_params$upper_bound)

#### Build timemap -------------------------------------------------------------

# pull timemap
tm <- con %>% tbl("tm") %>% collect()

# restrict time_map to distinct visit types & period before
tm <- tm %>%
  inner_join(select(index_dx_dates,patient_id,index_date),by = "patient_id") %>%    # filter to enrolled ids
  mutate(days_since_index = svcdate-index_date) %>%
  filter(between(days_since_index,-delay_params$upper_bound,0)) %>%
  distinct(patient_id,svcdate,days_since_index,outpatient,ed,inpatient,other)

save(tm, file = paste0(sim_out_path,"/delay_tm.RData"))
# load(paste0(sim_out_path,"/delay_tm.RData"))

#### Generate obs values for simulation ----------------------------------------
sim_obs <- tm %>% 
  filter(!(outpatient==0 & ed==0 & inpatient==0 & other == 0)) %>% 
  distinct(patient_id,days_since_index) %>% 
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>% 
  mutate(obs=row_number())

#### Generate Counts of Each ICD-9 & 10 ----------------------------------------

## Filter to final counts 
all_dx_visits <- con %>%
  tbl("dx_dates_alg1") %>%
  collect() %>%
  inner_join(select(index_dx_dates,patient_id,index_date),by = "patient_id") %>%    # filter to enrolled ids
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>% 
  select(patient_id,dx,dx_ver,inpatient,date,days_since_index)

# make sure to only count cases with continuous enrollment (i.e., join with filtered index dates)
all_dx_visits <- all_dx_visits %>%
  inner_join(distinct(index_dx_dates,patient_id), by = "patient_id")

# require corresponding obs from tm  (note: two steps so we keep the day 0 visits)
all_dx_visits <- bind_rows(all_dx_visits %>% 
            filter(days_since_index>=0),
          all_dx_visits %>% 
            inner_join(distinct(sim_obs,patient_id,days_since_index))) %>% 
  arrange(patient_id,days_since_index)

# dx counts
dx_counts <- all_dx_visits %>%
  distinct(patient_id,dx,dx_ver,days_since_index) %>%
  count(dx,dx_ver,days_since_index)

# visit counts
visit_counts <- all_dx_visits %>%
  distinct(patient_id,dx_ver,days_since_index) %>%
  count(dx_ver,days_since_index)

# populate missing valules in visit counts (i.e., assign 0 to days missing)
visit_counts <- tibble(days_since_index=-delay_params$upper_bound:delay_params$upper_bound) %>%
  mutate(dx_ver=map(days_since_index,~unique(dx_counts$dx_ver))) %>%
  unnest(dx_ver) %>%
  arrange(dx_ver,days_since_index) %>%
  left_join(visit_counts,by = c("days_since_index", "dx_ver")) %>%
  mutate(n = replace_na(n,0))

dx_counts <- dx_counts %>%
  rename(n_dx=n) %>%
  inner_join(visit_counts,by = c("dx_ver", "days_since_index")) %>%
  rename(n_visits = n) %>%
  mutate(frac = 100*n_dx/n_visits) # fraction of visits on a given day 

# save version for cluster analysis
save(dx_counts,visit_counts,file = paste0(cluster_out_path,"/dx_counts.RData"))
# load(paste0(cluster_out_path,"/dx_counts.RData"))

#### Save all dx visits for simulation -----------------------------------------

all_dx_visits <- all_dx_visits %>%
  filter(days_since_index<=0) %>%
  distinct(patient_id,dx,dx_ver,days_since_index)

tmp <- all_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index)

visit_counts <- bind_rows(tmp,visit_counts %>%
            filter(days_since_index<=0))

save(index_dx_dates,file = paste0(sim_out_path,"/index_dx_dates.RData"))
save(sim_obs,file = paste0(sim_out_path,"/sim_obs.RData"))
save(visit_counts,file = paste0(sim_out_path,"/visit_counts.RData"))
save(all_dx_visits,file = paste0(sim_out_path,"/dx_visits.RData"))
save(all_dx_visits,visit_counts,index_dx_dates,sim_obs,file = paste0(sim_out_path,"/all_dx_visits.RData"))
# load(paste0(sim_out_path,"/all_dx_visits.RData"))


########################
#### Pull demo data ####
########################
  
demographic_info <- haven::read_sas(paste0("/Shared/AML/kaiser_data/", 
                                           str_remove(cond_name, "_kaiser"), "/", 
                                           str_remove(cond_name, "_kaiser"), "_demographics_12sep23final.sas7bdat"))
demo1 <- demographic_info %>% 
  rename(patient_id = STUDYID,
         race = race1,
         sex = sex_admin) %>% 
  mutate(dobyr = lubridate::year(dob)) %>% 
  inner_join(select(index_dx_dates,patient_id,index_date),by = "patient_id") %>% 
  mutate(index_year = year(as_date(index_date))) 

## Save demographic info
save(demo1, file = paste0(sim_out_path,"/demo_data.RData"))

##########################################
#### Pull in More Granular Visit Info ####
##########################################

## STDPLAC ---------------------------------------------------------------------

stdplac_visits <- con %>% tbl("stdplac_visits") %>% collect()

tm_stdplac <- select(tm,patient_id,svcdate,days_since_index) %>% 
  inner_join(stdplac_visits,by = join_by(patient_id, svcdate))

save(tm_stdplac, file = paste0(sim_out_path,"/visit_info.RData"))

