
rm(list = ls())
library(tidyverse)
library(smallDB)


con <- DBI::dbConnect(RSQLite::SQLite(), "/Shared/AML/truven_extracts/small_dbs/dengue/dengue.db")
# con <- DBI::dbConnect(RSQLite::SQLite(), "~/Data/MarketScan/truven_extracts/small_dbs/dengue/dengue.db")

# pull timemap
tm <- con %>% tbl("tm") %>% collect()

##########################
#### Script Functions ####
##########################

# get inpatient core data
get_in_core <- function(source,year){
  con %>%
    tbl(glue::glue("inpatient_core_{source}_{year}")) %>% 
    select(any_of(c("caseid","admdate","disdate","los","patient_id"))) %>%
    collect()
}

# get inpatient diagnosis data
get_dx_data_in <- function(source,year){
  if (as.integer(year) > 14){
    tmp1 <- con %>%
      tbl(glue::glue("inpatient_core_{source}_{year}")) %>%
      select(any_of(c("caseid","admdate","disdate","los","patient_id"))) %>%
      inner_join(tbl(con,glue::glue("inpatient_dx9_{source}_{year}")),by="caseid") %>%
      distinct(patient_id,caseid,admdate,dx) %>%
      collect() %>%
      mutate(dx_ver=9L,
             patient_id=patient_id)
    
    tmp2 <- con %>%
      tbl(glue::glue("inpatient_core_{source}_{year}")) %>%
      select(any_of(c("caseid","admdate","disdate","los","patient_id"))) %>%
      inner_join(tbl(con,glue::glue("inpatient_dx10_{source}_{year}")),by="caseid") %>%
      distinct(patient_id,caseid,admdate,dx) %>%
      collect() %>%
      mutate(dx_ver=10L,
             patient_id=patient_id)
    
    rbind(tmp1,tmp2)
    
  } else {
    con %>%
      tbl(glue::glue("inpatient_core_{source}_{year}")) %>%
      select(any_of(c("caseid","admdate","disdate","los","patient_id"))) %>%
      inner_join(tbl(con,glue::glue("inpatient_dx_{source}_{year}")),by="caseid") %>%
      distinct(patient_id,caseid,admdate,dx) %>%
      collect() %>%
      mutate(dx_ver=9L,
             patient_id=patient_id)
  }
  
}

# get facility headers data
get_dx_data_fac <-function(source,year){
  
  tbl(con,glue::glue("facility_dx_{source}_{year}")) %>%
    inner_join(tbl(con,glue::glue("facility_core_{source}_{year}")) %>% 
                 select(fachdid,caseid), by = "fachdid") %>% 
    distinct(patient_id,svcdate,caseid,dx,dx_ver) %>%
    collect() 
  
}


##################################
#### Pull Inpatient Diagnoses ####
##################################

## Collect base inpatient data -------------------------------------------------

# core information
inpatient_core <- collect_plan(con) %>%
  mutate(data=map2(source,year,get_in_core)) %>%
  select(data) %>%
  unnest(cols = c(data)) %>% 
  mutate(disdate=ifelse(is.na(disdate),admdate+los,disdate))

# diagnoses
inpatient_dx <- collect_plan(con) %>%
  mutate(data=map2(source,year,get_dx_data_in)) %>%
  select(data) %>%
  unnest(cols = c(data))

## Collect inpatient services data ---------------------------------------------

in_serv_dx <- con %>% 
  tbl("inpatient_services_dx") %>% 
  select(patient_id,caseid,dx,dx_ver,svcdate,admdate,disdate) %>% 
  collect()

## Collect facility headers data -----------------------------------------------

facility_dx <- collect_plan(con) %>%
  filter(year!="01") %>% 
  mutate(data=map2(source,year,get_dx_data_fac)) %>% 
  mutate(data=map(data,~mutate(.,caseid = as.integer(caseid)))) %>% 
  select(data) %>%
  unnest(cols = c(data))

# partition into inpatient versus outpatient
in_facility_dx <- filter(facility_dx,!is.na(caseid) & 
                           caseid!=0)

out_facility_dx <- facility_dx %>% 
  anti_join(in_facility_dx,by = join_by(patient_id, svcdate, caseid, dx, dx_ver))


#############################################
#### Construct inpatient diagnosis dates ####
#############################################

## Motivation / Algorithm ------------------------------------------------------

### We use the following approach to constuct inpatient diagnosis dates for the
### delayed diagnosis simulation

# 1) Treat every inpatient date as a distinct observation. If we are able to 
#    assign a diagnosis date to a given diagnosis we try to do so.

# 2) Apply all dx to their respective date during the stay through the end of 
#    that stay...unless we can't assign a specific service date. Then apply to 
#    the whole stay

# To implement this procedure we apply the following logic
#  -If we have a specific service date in either the facility headers or 
#   the inpatient services table we use that...and apply foward until the stay 
#   ends.
#  -If we cannot find a specific service date for one of the diagnoses in a given
#   stay (i.e., from the main inpatient day). Then we apply that to the entire 
#   stay
#  -If we find a diagnoses service date in either the facility headers or 
#   inpatient services table that occurs before or after a given stay, we only 
#   apply to that particular date.

# Note: the current code doe not address situtions where a service date occurs 
#       before the official start date of a given stay for a condition that is 
#       is also diagnosed during the stay. For example suppose we have a service 
#       date for diagnoses X 2 days before a stay and then again on the second 
#       day of the stay. In this scenario we would end up having the diagnosis 
#       appear 2 days before, not on the day before or the first of the stay, 
#       and then again on the second day of the stay through the end of the stay.
#       However, one may argue that if you see it before and during the stay the
#       diagnosis should be carried forward.

## Pull out distinct inpatient diagnoses ---------------------------------------

# distinct inpatient diagnoses with admission and discharge dates
in_dx <- inpatient_dx %>% 
  inner_join(inpatient_core,by = join_by(patient_id, caseid, admdate))

## Split up facility diagnoses -------------------------------------------------

# facility inpatient that match on caseid
in_fac1 <- in_facility_dx %>% 
  inner_join(inpatient_core,by = join_by(patient_id, caseid)) %>% 
  filter(svcdate<=disdate & svcdate>=admdate)

# facility inpatient that match on time but not caseid
in_fac2 <-  in_facility_dx %>%
  anti_join(in_fac1,by = join_by(patient_id, svcdate, caseid, dx, dx_ver)) %>% 
  select(-caseid) %>% 
  inner_join(inpatient_core,by = join_by(patient_id),relationship = "many-to-many") %>% 
  filter(svcdate<=disdate & svcdate>=admdate)

# remaining inpatient facility dx that match on neither time nor caseid
in_fac3 <- in_facility_dx %>% 
  anti_join(in_fac1,by = join_by(patient_id, svcdate, caseid, dx, dx_ver)) %>% 
  anti_join(select(in_fac2,-caseid),by = join_by(patient_id, svcdate, dx, dx_ver)) 

## Split up inpatient services -------------------------------------------------

# inpatient services that match on caseid
in_serv1 <- in_serv_dx %>% 
  inner_join(inpatient_core, by = join_by(patient_id, caseid, admdate, disdate))

# inpatient services that match on caseid but have a missing discharge date
in_serv2 <- in_serv_dx %>% 
  anti_join(in_serv1,by = join_by(patient_id, caseid, dx, dx_ver, svcdate, admdate, disdate)) %>% 
  select(-disdate) %>% 
  inner_join(inpatient_core,by = join_by(patient_id, caseid, admdate)) 


# NOTE: NEED TO LOOK AT BIGGER DATASET FOR CASES BEYOND THIS TO SEE IF OCCUR
# inpatient services that do not match on caseid but fall during a given stay
in_serv3 <- in_serv_dx %>% 
  anti_join(in_serv1,join_by(patient_id, caseid, dx, dx_ver, svcdate, admdate, disdate)) %>% 
  select(-disdate) %>% 
  anti_join(in_serv2,by = join_by(patient_id, caseid, dx, dx_ver, svcdate, admdate)) %>% 
  select(-caseid,-admdate) %>% 
  inner_join(inpatient_core,by = join_by(patient_id)) %>% 
  filter(svcdate<=disdate & svcdate>=admdate)

# inpatient services that cannot be matched to a stay
in_serv4 <- in_serv_dx %>% 
  anti_join(in_serv1,join_by(patient_id, caseid, dx, dx_ver, svcdate, admdate, disdate)) %>% 
  select(-disdate) %>% 
  anti_join(in_serv2,by = join_by(patient_id, caseid, dx, dx_ver, svcdate, admdate)) %>% 
  select(-caseid,-admdate) %>% 
  anti_join(in_serv3,by = join_by(patient_id, dx, dx_ver, svcdate)) 
  

## Split up the inpatient diagnoses --------------------------------------------

#remove diagnoses from inpatient dx that are contained in other sets

# inpatient facility headers where we can match a caseid
in_fac_visit_matched <- bind_rows(in_fac1,in_fac2) %>% 
  distinct(patient_id,caseid,dx,dx_ver,admdate)

# inpatient services where we can match a caseid
in_serv_visit_matched <- bind_rows(in_serv1,in_serv2,in_serv3) %>% 
  distinct(patient_id,caseid,dx,dx_ver,admdate)

# now remove the inpatient diagnoses with facility header dates
in_dx1 <- inpatient_dx %>% 
  anti_join(in_fac_visit_matched,by = join_by(patient_id, caseid, admdate, dx, dx_ver)) 

# now remove the inpatient diagnoses with inpatient services dates
in_dx2 <- in_dx1 %>% 
  anti_join(in_serv_visit_matched,by = join_by(patient_id, caseid, admdate, dx, dx_ver))

# the remaing diagnoses that do not have dates in the facility header or inpatient services
in_dx_unmatched <- in_dx2
rm(in_dx1,in_dx2,in_fac_visit_matched,in_serv_visit_matched)
gc()

## Create the inpatient diagnosis spans ----------------------------------------

# remaining inpatient visit span
in_dx_span <- in_dx_unmatched %>% 
  inner_join(inpatient_core, by = join_by(patient_id, caseid, admdate)) %>% 
  distinct(patient_id,dx,dx_ver,start = admdate,end=disdate) %>% 
  mutate(date = map2(start,end,~.x:.y)) %>% 
  unnest(date) %>% 
  distinct(patient_id,dx,dx_ver,date)

# facility header span
in_fac_span <- bind_rows(in_fac1,in_fac2) %>% 
  distinct(patient_id,svcdate,caseid,dx,dx_ver,admdate) %>% 
  inner_join(inpatient_core, by = join_by(patient_id, caseid, admdate)) %>% 
  distinct(patient_id,dx,dx_ver,start = svcdate,end=disdate) %>% 
  mutate(date = map2(start,end,~.x:.y)) %>% 
  unnest(date) %>% 
  distinct(patient_id,dx,dx_ver,date)

# inpatient services span
in_serv_span <- bind_rows(in_serv1,in_serv2,in_serv3) %>% 
  distinct(patient_id,svcdate,caseid,dx,dx_ver,admdate) %>% 
  inner_join(inpatient_core, by = join_by(patient_id, caseid, admdate)) %>% 
  distinct(patient_id,dx,dx_ver,start = svcdate,end=disdate) %>% 
  mutate(date = map2(start,end,~.x:.y)) %>% 
  unnest(date) %>% 
  distinct(patient_id,dx,dx_ver,date)


# unspanned diagnoses...do no match to stay
bind_rows(in_fac3,in_serv4)

## Combine diagnoses
in_dx_dates <- bind_rows(in_dx_span,
          in_fac_span,
          in_serv_span,
          distinct(bind_rows(in_fac3,in_serv4),patient_id,dx,dx_ver,date=svcdate)) %>% 
  distinct()


# check that we still have all the same diagnoses across patients
bind_rows(distinct(in_facility_dx,patient_id,dx,dx_ver),
          distinct(inpatient_dx,patient_id,dx,dx_ver),
          distinct(in_serv_dx,patient_id,dx,dx_ver)) %>% 
  distinct() %>% 
  nrow()
  
in_dx_dates %>% 
  distinct(patient_id,dx,dx_ver) %>% 
  nrow()

## Cleanup ---------------------------------------------------------------------
rm(in_dx,
   in_fac1,in_fac2,in_fac3,
   in_serv1,in_serv2,in_serv3,in_serv4,
   in_dx_unmatched,in_dx_span,in_fac_span,in_serv_span)
gc()




##################################################
#### Construct the outpatient diagnosis dates ####
##################################################