

library(tidyverse)
library(icd)
library(codeBuildr)


cond_name <- "sepsis_pre_covid"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)

sim_res_path <- paste0(delay_params$out_path,"sim_results/exponential_cp14/")

###################
#### Load Data ####
###################

## Load index cases ------------------------------------------------------------
load(paste0(delay_params$out_path,"index_cases.RData"))
# index_cases

## Load in prebuilt demographic data -------------------------------------------
load(paste0(delay_params$base_path,"delay_results/demo_data.RData"))
# demo1
# demo2
# rural_visits

demo1 <- demo1 %>% inner_join(select(index_cases,patient_id), by = "patient_id")
demo2 <- demo2 %>% inner_join(select(index_cases,patient_id), by = "patient_id")
rural_visits <- rural_visits %>% inner_join(select(index_cases,patient_id), by = "patient_id")

## load in dx visits -----------------------------------------------------------
load(paste0(delay_params$base_path,"delay_results/all_dx_visits.RData"))
# visit_counts
# all_dx_visits
all_dx_visits <- all_dx_visits %>% inner_join(select(index_cases,patient_id), by = "patient_id")
sim_obs <- sim_obs %>% inner_join(select(index_cases,patient_id), by = "patient_id")
rm(visit_counts)

## load in timemap -------------------------------------------------------------
load(paste0(delay_params$base_path,"delay_results/delay_tm.RData"))
# tm
tm <- tm %>% inner_join(select(index_cases,patient_id), by = "patient_id")

## load simulation results -----------------------------------------------------
load(paste0(sim_res_path,"sim_res_ssd.RData"))
# sim_res_ssd

## reduce sim_obs --------------------------------------------------------------
tmp <- sim_res_ssd %>% 
  mutate(obs = map(res,~distinct(.,obs))) %>% 
  select(obs) %>% 
  unnest(obs) %>% 
  distinct(obs)

sim_obs <- inner_join(tmp,sim_obs) 

## load prescriptions ----------------------------------------------------------
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path,"sepsis_revised10.db"))

rx_visits <- con %>% tbl("all_rx_visits") %>% collect()

rx_visits <- rx_visits %>% 
  inner_join(select(index_cases,patient_id,index_date)) %>% 
  filter(date<=index_date & date>=(index_date-180))

gc()

## Load Bootstrap sample -------------------------------------------------------
load(paste0(delay_params$out_path,"sim_results/boot_data.RData"))
boot_data$boot_sample

############################
#### Prepare Visit Info ####
############################

index_locations <- tm %>% 
  filter(days_since_index==0) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,setting_type,admdate) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  mutate(dow=weekdays(as_date(admdate))) %>% 
  mutate(weekend = dow %in% c("Saturday","Sunday")) %>% 
  select(patient_id,setting,weekend,dow) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(tribble(~outpatient,~ed,~inpatient,~setting_label,
                     T,F,F,"Out only",
                     T,T,F,"Out and ED",
                     T,T,T,"All three",
                     T,F,T,"Out and inpatient",
                     F,T,F,"ED only",
                     F,T,T,"ED and inpatient",
                     F,F,T,"Inpatient only",
                     F,F,F,"none")) %>% 
  mutate(setting_label = fct_relevel(setting_label,
                                     "Out only",
                                     "ED only",
                                     "Inpatient only",
                                     "Out and ED",
                                     "Out and inpatient",
                                     "ED and inpatient",
                                     "All three",
                                     "none")) 

index_locations %>% count(dow,ed)

obs_locations <- tm %>% 
  left_join(sim_obs,by = c("patient_id", "days_since_index")) %>% 
  filter(!is.na(obs)) %>% 
  filter(setting_type !=4) %>%
  distinct(patient_id,obs,setting_type) %>% 
  mutate(setting=smallDB::setting_type_labels(setting_type)) %>%
  select(patient_id,obs,setting) %>%
  mutate(value = TRUE) %>% 
  spread(key= setting, value = value,fill = FALSE) %>% 
  inner_join(tribble(~outpatient,~ed,~inpatient,~setting_label,
                     T,F,F,"Out only",
                     T,T,F,"Out and ED",
                     T,T,T,"All three",
                     T,F,T,"Out and inpatient",
                     F,T,F,"ED only",
                     F,T,T,"ED and inpatient",
                     F,F,T,"Inpatient only",
                     F,F,F,"none")) %>% 
  mutate(setting_label = fct_relevel(setting_label,
                                     "Out only",
                                     "ED only",
                                     "Inpatient only",
                                     "Out and ED",
                                     "Out and inpatient",
                                     "ED and inpatient",
                                     "All three",
                                     "none")) 

reg_demo <- demo1 %>% 
  mutate(female=(sex==2),
         age = index_year-dobyr,
         stdrace = as.numeric(stdrace)) %>% 
  left_join(tibble(stdrace = c(0,1,2,4,9),
                   race = c("Missing/Unknown","White","Black","Hispanic","Other")),
            by = "stdrace") %>% 
  mutate(race = fct_relevel(race,"White"))

reg_demo <- reg_demo %>% 
  left_join(demo2 %>% 
              filter(index_date<=dtend & index_date>=dtstart) %>% 
              mutate(msa_new = msa %in% c("0","")) %>% 
              mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
              mutate(source = as.factor(source)) %>% 
              distinct(patient_id,source,msa=msa_new),
            by = "patient_id")

age_cats <- c(-1,17,34,44,54,64,130)

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))


### Prep weekend visit info ------------

# # add weekend and demo to sim data
# sim_res <- sim_res %>% 
#   select(obs,trial) %>% 
#   inner_join(obs_locations) %>% 
#   inner_join(sim_res_sim_obs %>% 
#                inner_join(select(index_dx_dates,patient_id,index_date), by = "patient_id") %>% 
#                mutate(vis_date = index_date+days_since_index) %>% 
#                mutate(dow=weekdays(as_date(vis_date))) %>% 
#                select(obs,dow) %>% 
#                mutate(weekend = dow %in% c("Saturday","Sunday")), by = "obs")


### Prepare Drug Indicators ----------------------------------------------------
red_book <- codeBuildr::load_redbook()

tmp1 <- red_book %>% filter(str_detect(tolower(THRGRDS),"immuno")) %>% .$NDCNUM
tmp2 <- red_book %>% filter(str_detect(tolower(THRDTDS),"prednisone")) %>% .$NDCNUM

immuno_rx_ids <- rx_visits %>% 
  filter(ndcnum %in% tmp1) %>% 
  filter(date<index_date-14 & (date+daysupp)>=(index_date-30)) %>% 
  distinct(patient_id) %>% 
  mutate(immunosuppressant=1L)


pred_rx_ids <- rx_visits %>% 
  filter(ndcnum %in% tmp2) %>% 
  filter(date<index_date-14 & (date+daysupp)>=(index_date-30)) %>% 
  distinct(patient_id) %>% 
  mutate(prednisone=1L)

reg_demo <- reg_demo %>% 
  left_join(immuno_rx_ids) %>% 
  left_join(pred_rx_ids) %>% 
  mutate_at(vars(immunosuppressant,prednisone),~replace_na(.,0L))


## Prepare comorbidity indicators ----------------------------------------------
tmp1 <- all_dx_visits %>% 
  filter(days_since_index < -14) %>% 
  filter(dx_ver==9) %>% 
  select(patient_id,dx)

tmp2 <- all_dx_visits %>% 
  filter(days_since_index < -14) %>% 
  filter(dx_ver==10) %>% 
  select(patient_id,dx)

tmp1 <- icd9_comorbid_elix(tmp1, 
                           visit_name = "patient_id",
                           icd_name = "dx",
                           return_df = TRUE) %>% 
  as_tibble() %>% 
  gather(key = key, value = value, -patient_id) %>% 
  filter(value == TRUE)

tmp2 <- icd10_comorbid_elix(tmp2, 
                           visit_name = "patient_id",
                           icd_name = "dx",
                           return_df = TRUE) %>% 
  as_tibble() %>% 
  gather(key = key, value = value, -patient_id) %>% 
  filter(value == TRUE)

tmp <- bind_rows(tmp1,tmp2) %>% 
  distinct() 

cm_count <- tmp %>% group_by(patient_id) %>% summarise(cm_count=n()) %>% mutate(patient_id = as.integer(patient_id))

cm_inds <- tmp %>% spread(key = key, value = value)  %>% mutate(patient_id = as.integer(patient_id))

reg_demo <- reg_demo %>% 
  left_join(cm_count) %>% 
  left_join(cm_inds) %>% 
  mutate_at(vars(cm_count:WeightLoss), ~replace_na(.,0L))


## Add Mental Health Codes -----------------------------------------------------
tmp1 <- all_dx_visits %>% 
  filter(days_since_index < -14) %>% 
  filter(dx_ver==9) %>% 
  select(patient_id,dx)

tmp2 <- all_dx_visits %>% 
  filter(days_since_index < -14) %>% 
  filter(dx_ver==10) %>% 
  select(patient_id,dx)


tmp1 <- tmp1 %>% 
  inner_join(codeBuildr::ccs9_mappings %>% 
               filter(ccs_code %in% c(651,653,654,657:661)) %>% 
               select(dx = icd9cm,ccs_code)) %>% 
  distinct(patient_id,ccs_code)

tmp2 <- tmp2 %>% 
  inner_join(codeBuildr::ccs10_mappings %>% 
               filter(ccs_code %in% c(651,653,654,657:661)) %>% 
               select(dx = icd10cm,ccs_code)) %>% 
  distinct(patient_id,ccs_code)

mh_inds <- bind_rows(tmp1,tmp2) %>% 
  distinct() %>% 
  mutate(ccs_code = paste0("ccs_",ccs_code)) %>% 
  mutate(value = TRUE) %>% 
  spread(key = ccs_code, value = value)


reg_demo <- reg_demo %>% 
  left_join(mh_inds) %>% 
  mutate_at(vars(ccs_651:ccs_661), ~replace_na(.,0L))


reg_demo %>% glimpse()

##############################
#### Save regression data ####
##############################

save(reg_demo,obs_locations,index_locations, file = paste0(delay_params$out_path,"sim_results/reg_vars.RData"))

reg_demo %>% glimpse()

################################
#### Run Miss Patient Model ####
################################
load(paste0(sim_res_path,"sim_res_ssd.RData"))
load(paste0(delay_params$out_path,"sim_results/reg_vars.RData"))
load(paste0(delay_params$base_path,"delay_results/all_dx_visits.RData"))
tmp <- sim_res_ssd %>%
  mutate(obs = map(res,~distinct(.,obs))) %>%
  select(obs) %>%
  unnest(obs) %>%
  distinct(obs)
sim_obs <- inner_join(tmp,sim_obs)
load(paste0(delay_params$out_path,"sim_results/boot_data.RData"))

rm(list = ls()[!(ls() %in% c("reg_demo","obs_locations","index_locations",
                             "boot_data","sim_res_ssd","delay_params","sim_obs",
                             "sim_res_path"))])

## Run Model 1 -----------------------------------------------------------------
gc()

tmp <- sim_res_ssd %>% 
  unnest(res) %>% 
  inner_join(sim_obs) 

tmp <- tmp %>% 
  group_by(sim_trial,boot_trial) %>% 
  nest() %>% 
  mutate(data = map(data,~distinct(.,patient_id,boot_id)))

sim_res <- tmp

rm(tmp,sim_res_ssd)
gc()

boot_data <- boot_data %>% select(boot_trial,boot_sample) %>% arrange(boot_trial)

run_model <- function(boot_trial,sim_data){
  
  reg_data <- boot_data$boot_sample[[boot_trial]] %>% 
    left_join(sim_data %>% 
                mutate(miss=1L),
              join_by(patient_id, boot_id)) %>%
    mutate(miss = replace_na(miss,0L)) %>% 
    left_join(reg_demo,join_by(patient_id))
  
  fit <- glm(miss~., 
             data = select(reg_data,miss,female,age_cat,source,month,immunosuppressant,prednisone,Alcohol:ccs_661),
             family = "binomial")
  
  bind_cols(broom::tidy(fit) %>% 
              mutate(or = exp(estimate)) %>% 
              select(term,or),
            exp(broom::confint_tidy(fit,func = stats::confint.default)))
  
}

cluster <- parallel::makeCluster(30)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("run_model","boot_data","reg_demo"),
                        envir=environment())


sim_res <- map(1:nrow(sim_res),~list(boot_trial = sim_res$boot_trial[[.]], data = sim_res$data[[.]]))

model_res <- parallel::parLapply(cl = cluster,
                                 sim_res,
                                 function(x){run_model(boot_trial = x$boot_trial,
                                                       sim_data = x$data)})
bind_rows(model_res) %>% 
  group_by(term) %>% 
  summarise_at(vars(or:conf.high),mean) %>% 
  write_csv(paste0(sim_res_path,"model1_estimate.csv"))

### Run model 2 ----------------------------------------------------------------

run_model2 <- function(boot_trial,sim_data){
  
  reg_data <- boot_data$boot_sample[[boot_trial]] %>% 
    left_join(sim_data %>% 
                mutate(miss=1L),
              join_by(patient_id, boot_id)) %>%
    mutate(miss = replace_na(miss,0L)) %>% 
    left_join(reg_demo,join_by(patient_id))
  
  fit <- glm(miss~., 
             data = select(reg_data,miss,female,age_cat,source,month,immunosuppressant,prednisone,cm_count,ccs_651:ccs_661),
             family = "binomial")
  
  bind_cols(broom::tidy(fit) %>% 
              mutate(or = exp(estimate)) %>% 
              select(term,or),
            exp(broom::confint_tidy(fit,func = stats::confint.default)))
  
}

parallel::clusterExport(cluster,c("run_model2"),envir=environment())

model_res2 <- parallel::parLapply(cl = cluster,
                                 sim_res,
                                 function(x){run_model2(boot_trial = x$boot_trial,
                                                       sim_data = x$data)})

bind_rows(model_res2) %>% 
  group_by(term) %>% 
  summarise_at(vars(or:conf.high),mean) %>% 
  write_csv(paste0(sim_res_path,"model2_estimate.csv"))


## Export Results --------------------------------------------------------------

mod1 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model1_estimate.csv")
mod2 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model2_estimate.csv")

mod1 %>% 
  mutate_at(vars(or:conf.high),~round(.,2))
  mutate(model1 = paste0())
  
  
  


#################################
#### Run Miss Duration Model ####
#################################
rm(list = ls()[!(ls() %in% c("reg_demo","obs_locations","index_locations",
                             "boot_data","sim_res_ssd","delay_params","sim_obs",
                             "sim_res_path"))])

sim_res <- sim_res_ssd %>% 
  mutate(res = map(res,~inner_join(.,sim_obs,by = join_by(obs)) %>% 
                     group_by(patient_id,boot_id) %>% 
                     summarise(duration = -min(days_since_index),.groups = "drop")))


run_model <- function(boot_trial,sim_data){
  
  reg_data <- boot_data$boot_sample[[boot_trial]] %>% 
    left_join(sim_data, join_by(patient_id, boot_id)) %>%
    mutate(duration = replace_na(duration,0L)) %>% 
    left_join(reg_demo,join_by(patient_id))
  
  fit <- lm(duration~., data = select(reg_data,duration,female,age_cat,source,month,immunosuppressant,prednisone,Alcohol:ccs_661))
  
  bind_cols(broom::tidy(fit) %>% 
              select(term,estimate),
            broom::confint_tidy(fit,func = stats::confint.default))
  
}

sim_res <- map(1:nrow(sim_res),~list(boot_trial = sim_res$boot_trial[[.]], data = sim_res$res[[.]]))

cluster <- parallel::makeCluster(30)

parallel::clusterCall(cluster, function() library(tidyverse))
parallel::clusterExport(cluster,c("run_model","boot_data","reg_demo"),
                        envir=environment())


model_res <- parallel::parLapply(cl = cluster,
                                 sim_res,
                                 function(x){run_model(boot_trial = x$boot_trial,
                                                       sim_data = x$data)})
bind_rows(model_res) %>% 
  group_by(term) %>% 
  summarise_at(vars(estimate:conf.high),mean) %>% 
  write_csv(paste0(sim_res_path,"model1_duration_estimate.csv"))

## Run second duration model ---------------------------------------------------

run_model2 <- function(boot_trial,sim_data){
  
  reg_data <- boot_data$boot_sample[[boot_trial]] %>% 
    left_join(sim_data, join_by(patient_id, boot_id)) %>%
    mutate(duration = replace_na(duration,0L)) %>% 
    left_join(reg_demo,join_by(patient_id))
  
  fit <- lm(duration~., data = select(reg_data,duration,female,age_cat,source,month,immunosuppressant,prednisone,cm_count,ccs_651:ccs_661))
  
  bind_cols(broom::tidy(fit) %>% 
              select(term,estimate),
            broom::confint_tidy(fit,func = stats::confint.default))
  
}

parallel::clusterExport(cluster,c("run_model2"),envir=environment())

model_res2 <- parallel::parLapply(cl = cluster,
                                  sim_res,
                                  function(x){run_model2(boot_trial = x$boot_trial,
                                                         sim_data = x$data)})

bind_rows(model_res2) %>% 
  group_by(term) %>% 
  summarise_at(vars(estimate:conf.high),mean) %>% 
  write_csv(paste0(sim_res_path,"model2_duration_estimate.csv"))



##############################
#### Run Miss Visit Model ####
##############################

library(tidyverse)
dur1 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model1_duration_estimate.csv")
dur2 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model2_duration_estimate.csv")
miss_pat1 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model1_estimate.csv")
miss_pat2 <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model2_estimate.csv")


rmarkdown::render("projects/sepsis_revised10/pre_covid/risk_report.Rmd", 
                  output_file = "/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/risk_report.html")
