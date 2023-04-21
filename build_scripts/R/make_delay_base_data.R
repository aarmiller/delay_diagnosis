

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
# cond_name <- "tb"

# Load Delay Params
load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name)

# connect to database
con <- DBI::dbConnect(RSQLite::SQLite(), paste0(delay_params$path,cond_name,".db"))

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

# Filter to enrolids that can be observed 
index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

index_dx_dates <- index_dx_dates %>%
  filter(time_before_index>=delay_params$upper_bound)

#### Build timemap -------------------------------------------------------------

# find available years and medicaid years
avail_data <- tibble(table = DBI::dbListTables(con)) %>%
  filter(str_detect(table,"inpatient_core")) %>%
  mutate(tmp = str_remove(table,"inpatient_core_")) %>%
  mutate(year = as.integer(str_sub(tmp,start = -2, end = -1))) %>%
  mutate(source = str_remove(tmp,"_[0-9]*")) %>%
  select(year,source)

medicaid_years <- filter(avail_data,source == "medicaid") %>%
  distinct(year) %>%
  arrange(year) %>%
  .$year

years <- filter(avail_data,source == "ccae") %>%
  distinct(year) %>%
  arrange(year) %>%
  .$year

collect_tab <- collect_table(years = years,medicaid_years = medicaid_years)

# pull timemap
tm <- build_time_map(db_con = con,
                     collect_tab = collect_tab)


# restrict time_map to distinct visit types & period before
tm <- tm %>%
  inner_join(select(index_dx_dates,enrolid,index_date),by = "enrolid") %>%    # filter to enrolled ids
  mutate(days_since_index = admdate-index_date) %>%
  filter(between(days_since_index,-delay_params$upper_bound,0)) %>%
  distinct(enrolid,admdate,disdate,days_since_index,stdplac,setting_type)

save(tm, file = paste0(sim_out_path,"/delay_tm.RData"))
# load(paste0(sim_out_path,"/delay_tm.RData"))

#### Generate obs values for simulation ----------------------------------------
sim_obs <- tm %>% 
  distinct(enrolid,days_since_index) %>% 
  filter(between(days_since_index,-delay_params$upper_bound,-1)) %>% 
  mutate(obs=row_number())

#### Generate Counts of Each ICD-9 & 10 ----------------------------------------

## Filter to final counts 
all_dx_visits <- con %>%
  tbl("all_dx_visits") %>%
  collect() %>%
  filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound))

# make sure to only count cases with continuous enrollment (i.e., join with filtered index dates)
all_dx_visits <- all_dx_visits %>%
  inner_join(distinct(index_dx_dates,enrolid), by = "enrolid")

# dx counts
dx_counts <- all_dx_visits %>%
  distinct(enrolid,dx,dx_ver,days_since_index) %>%
  count(dx,dx_ver,days_since_index)

# visit counts
visit_counts <- all_dx_visits %>%
  distinct(enrolid,dx_ver,days_since_index) %>%
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
  distinct(enrolid,dx,dx_ver,days_since_index)

tmp <- all_dx_visits %>%
  distinct(enrolid,days_since_index) %>%
  count(days_since_index)

visit_counts <- bind_rows(tmp,visit_counts %>%
            filter(days_since_index<=0))


save(all_dx_visits,visit_counts,index_dx_dates,sim_obs,file = paste0(sim_out_path,"/all_dx_visits.RData"))
# load(paste0(sim_out_path,"/all_dx_visits.RData"))


########################
#### Pull demo data ####
########################

### Get Basic Demographics -----------------------------------------------------

med_collect <- collect_tab %>% filter(source=="medicaid") %>% 
  distinct(source,year)

ccae_mdcr_collect <- collect_tab %>% filter(source!="medicaid") %>% 
  distinct(source,year)

# collect medicaid
tmp_med <- smallDB::gether_table_data(collect_tab = med_collect,
                          table = "enrollees",
                          vars = c("enrolid","sex","dobyr","stdrace"),
                          db_con = con,
                          collect_n = Inf) %>% 
  select(year,table_data) %>% 
  mutate(year = as.integer(year)+2000L) %>% 
  unnest(table_data) %>% 
  inner_join(select(index_dx_dates,enrolid,index_date)) %>% 
  mutate(index_year = year(as_date(index_date))) %>% 
  filter(index_year==year) 

# collect ccae medicare
tmp_ccae_mdcr <- smallDB::gether_table_data(collect_tab = ccae_mdcr_collect,
                                            table = "enrollees",
                                            vars = c("enrolid","sex","dobyr"),
                                            db_con = con,
                                            collect_n = Inf) %>% 
  select(year,table_data) %>% 
  mutate(year = as.integer(year)+2000L) %>% 
  unnest(table_data) %>% 
  inner_join(select(index_dx_dates,enrolid,index_date)) %>% 
  mutate(index_year = year(as_date(index_date))) %>% 
  filter(index_year==year) %>% 
  distinct()

demo1 <- bind_rows(tmp_ccae_mdcr,tmp_med)
rm(tmp_ccae_mdcr,tmp_med)

# collect state, msa and egeoloc by enroll_period
tmp <- smallDB::gether_table_data(collect_tab = collect_tab,
                                    table = "enrollment_detail",
                                    vars = c("enrolid","dtend","dtstart","msa","egeoloc"),
                                    db_con = con,collect_n = Inf)

tmp <- tmp %>% 
  select(source,year,table_data) %>% 
  mutate(table_data = map(table_data,~mutate_all(.,as.character))) %>% 
  unnest(table_data) %>% 
  mutate_at(vars(enrolid,dtend,dtstart,egeoloc),as.integer)

# filter to delay observation window
demo2 <- tmp %>% 
  inner_join(select(index_dx_dates,enrolid,index_date),by = "enrolid") %>% 
  filter(dtstart<=index_date & dtend>=(index_date-delay_params$upper_bound)) %>%
  select(-year)


## Load in rurality info -------------------------------------------------------

rural_visits <- collect_tab %>% 
  distinct(source,year) %>% 
  mutate(data = map2(source,year,~tbl(con,paste0("outpatient_core_",.x,"_",.y)) %>% 
                       filter(stdplac==72) %>% 
                       distinct(enrolid,svcdate) %>% 
                       collect()))

rural_visits <- rural_visits %>% 
  unnest(data) %>% 
  distinct(enrolid,svcdate)


## Save demographic info
save(demo1,demo2,rural_visits, file = paste0(sim_out_path,"/demo_data.RData"))


