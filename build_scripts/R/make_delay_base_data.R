
rm(list = ls())
library(tidyverse)
library(bit64)
# library(parallel)
library(smallDB)
# library(delaySim)
library(lubridate)

############################################
#### Load Data and Prepare Output paths ####
############################################

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "pd"

# Load Delay Params
load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")

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
  distinct(patient_id,svcdate,days_since_index,mdcr,ccae,medicaid,outpatient,ed,obs_stay,inpatient,rx)

save(tm, file = paste0(sim_out_path,"/delay_tm.RData"))
# load(paste0(sim_out_path,"/delay_tm.RData"))

#### Generate obs values for simulation ----------------------------------------
sim_obs <- tm %>% 
  filter(!(outpatient==0 & ed==0 & obs_stay==0 & inpatient==0)) %>% 
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
  mutate(dx_ver=map(days_since_index,~c(9,10))) %>%
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

### Get Basic Demographics -----------------------------------------------------

med_collect <- collect_plan(con) %>% filter(source=="medicaid") %>% 
  distinct(source,year)

ccae_mdcr_collect <- collect_plan(con) %>% filter(source!="medicaid") %>% 
  distinct(source,year)

# collect medicaid
if (nrow(med_collect)>0){
  tmp_med <- smallDB::gether_table_data(collect_tab = med_collect,
                                        table = "enrollees",
                                        vars = c("patient_id","sex","dobyr","stdrace"),
                                        db_con = con,
                                        collect_n = Inf) %>% 
    select(year,table_data) %>% 
    mutate(year = as.integer(year)+2000L) %>% 
    mutate(table_data = map(table_data,~mutate(.,stdrace=as.character(stdrace)))) %>% 
    unnest(table_data) %>% 
    inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
    mutate(index_year = year(as_date(index_date))) %>% 
    filter(index_year==year) 
} else {
  tmp_med <- tibble()
}


# collect ccae medicare
tmp_ccae_mdcr <- smallDB::gether_table_data(collect_tab = ccae_mdcr_collect,
                                            table = "enrollees",
                                            vars = c("patient_id","sex","dobyr"),
                                            db_con = con,
                                            collect_n = Inf) %>% 
  select(year,table_data) %>% 
  mutate(year = as.integer(year)+2000L) %>% 
  unnest(table_data) %>% 
  inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
  mutate(index_year = year(as_date(index_date))) %>% 
  filter(index_year==year) %>% 
  distinct()

demo1 <- bind_rows(tmp_ccae_mdcr,tmp_med)
rm(tmp_ccae_mdcr,tmp_med)

# collect state, msa and egeoloc by enroll_period
tmp <- smallDB::gether_table_data(collect_tab = collect_plan(con),
                                  table = "enrollment_detail",
                                  vars = c("patient_id","dtend","dtstart","msa","egeoloc"),
                                  db_con = con,collect_n = Inf)

tmp <- tmp %>% 
  select(source,year,table_data) %>% 
  mutate(table_data = map(table_data,~mutate_all(.,as.character))) %>% 
  unnest(table_data) %>% 
  mutate_at(vars(patient_id,dtend,dtstart,egeoloc),as.integer)

# filter to delay observation window
demo2 <- tmp %>% 
  inner_join(select(index_dx_dates,patient_id,index_date),by = "patient_id") %>% 
  filter(dtstart<=index_date & dtend>=(index_date-delay_params$upper_bound)) %>%
  select(-year)


## Load in rurality info -------------------------------------------------------

rural_visits <- collect_plan(con) %>% 
  distinct(source,year) %>% 
  mutate(data = map2(source,year,~tbl(con,paste0("outpatient_core_",.x,"_",.y)) %>% 
                       filter(stdplac==72) %>% 
                       distinct(patient_id,svcdate) %>% 
                       collect()))

rural_visits <- rural_visits %>% 
  unnest(data) %>% 
  distinct(patient_id,svcdate)


## Save demographic info
save(demo1,demo2,rural_visits, file = paste0(sim_out_path,"/demo_data.RData"))


##########################################
#### Pull in More Granular Visit Info ####
##########################################

## STDPLAC ---------------------------------------------------------------------

stdplac_visits <- con %>% tbl("stdplac_visits") %>% collect()

tm_stdplac <- select(tm,patient_id,svcdate,days_since_index) %>% 
  inner_join(stdplac_visits,by = join_by(patient_id, svcdate))


## STDPROV ---------------------------------------------------------------------

stdprov_visits <- con %>% tbl("stdprov_visits") %>% collect()

tm_stdprov <- select(tm,patient_id,svcdate,days_since_index) %>% 
  inner_join(stdprov_visits,by = join_by(patient_id, svcdate))


## SVCSCAT ---------------------------------------------------------------------

svcscat_visits <- con %>% tbl("svcscat_visits") %>% collect()

tm_svcscat <- select(tm,patient_id,svcdate,days_since_index) %>% 
  inner_join(svcscat_visits,by = join_by(patient_id, svcdate))


save(tm_stdplac,tm_stdprov,tm_svcscat, file = paste0(sim_out_path,"/visit_info.RData"))


######################
#### Pull Caseids ####
######################

caseids <- con %>% 
  tbl("tm_full") %>% 
  filter(!is.na(caseid)) %>% 
  select(patient_id,svcdate,caseid,admdate,disdate) %>% 
  collect() 

caseids <- caseids %>% 
  distinct(patient_id,date=svcdate,caseid) %>% 
  inner_join(select(index_dx_dates,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-delay_params$upper_bound,0)) %>% 
  left_join(sim_obs,by = join_by(patient_id, days_since_index)) %>% 
  distinct(patient_id,date,caseid,days_since_index,obs)

save(caseids, file = paste0(sim_out_path,"/caseids.RData"))  
  