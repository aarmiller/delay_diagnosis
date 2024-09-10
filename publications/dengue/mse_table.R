
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)
library(changepoint)
library(trend)

params <- list(proj = "dengue")

# load delay_parms
load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[params$proj]]

project_test <- codeBuildr::avail_ssd_codes() %>% 
  filter(name == params$proj) %>% nrow()

if(project_test>0){
  
  cond_label <- codeBuildr::avail_ssd_codes() %>% 
    filter(name == params$proj) %>% 
    .$description
  
  ssd_codes <- codeBuildr::load_ssd_codes(params$proj) %>% 
    filter(type %in% c("icd9","icd10")) %>% 
    mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
    select(dx = code,dx_ver)
  
} else {
  cond_label <- stringr::str_to_title(stringr::str_replace(params$proj, "_", " "))
  
  ssd_codes <- codeBuildr::load_ssd_codes(stringr::str_split(params$proj, "_")[[1]][1]) %>% 
    filter(type %in% c("icd9","icd10")) %>% 
    mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
    select(dx = code,dx_ver)
}

delay_base_path <- paste0(delay_params$base_path,"delay_results/")

# identify test dates
load(paste0(delay_params$out_path,"index_cases.RData"))
load(paste0(delay_base_path,"all_dx_visits.RData"))
rm(visit_counts, sim_obs, index_dx_dates)

# update all_dx_visits
all_dx_visits <- all_dx_visits %>%
  inner_join(index_cases %>% select(patient_id, shift), by = "patient_id") %>% 
  mutate(days_since_index = days_since_index - shift) %>% 
  select(-shift) %>% 
  filter(days_since_index<=0) 

patient_ids <- index_cases %>% distinct(patient_id)

# visit counts
visit_counts <- all_dx_visits %>%
  distinct(patient_id,dx_ver,days_since_index) %>%
  count(dx_ver,days_since_index)

# populate missing values in visit counts (i.e., assign 0 to days missing)
visit_counts <- tibble(days_since_index=-delay_params$upper_bound:delay_params$upper_bound) %>%
  mutate(dx_ver=map(days_since_index,~c(9,10))) %>%
  unnest(dx_ver) %>%
  arrange(dx_ver,days_since_index) %>%
  left_join(visit_counts,by = c("days_since_index", "dx_ver")) %>%
  mutate(n = replace_na(n,0))

tmp <- all_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index)

visit_counts <- bind_rows(tmp,visit_counts %>%
                            filter(days_since_index<=0))


fit_model <- function(data,cp,model="lm",periodicity = FALSE, return_fit=FALSE){
  
  tmp <- data %>% 
    mutate(before_cp = period>=cp)
  
  if (model=="lm"){
    if (periodicity){
      fit <- lm(n~period*before_cp+dow, data = tmp)
    } else {
      fit <- lm(n~period*before_cp, data = tmp) 
    }
  } else if (model=="quad"){
    if (periodicity){
      fit <- lm(n~poly(period,2)*before_cp+dow, data = tmp)
    } else {
      fit <- lm(n~poly(period,2)*before_cp, data = tmp) 
    }
  } else if (model=="cubic") {
    if (periodicity){
      fit <- lm(n~poly(period,3)*before_cp+dow, data = tmp)
    } else {
      fit <- lm(n~poly(period,3)*before_cp, data = tmp) 
    }
  } else if (model=="exp"){
    if (periodicity){
      fit <- lm(n~log(period)*before_cp+dow, data = tmp)
    } else {
      fit <- lm(n~log(period)*before_cp, data = tmp) 
    }
  }
  
  rmse <- sqrt(mean(fit$residuals^2))
  
  if (return_fit){
    return(list(fit=fit,
                rmse=rmse))
  } else {
    return(rmse)
  }
  
}

find_cp <- function(data,cp_range,model="lm",periodicity = FALSE){
  out1 <- tibble(cp = cp_range) %>% 
    mutate(rmse = map_dbl(cp,~fit_model(data = data,
                                        cp = .,
                                        model = model,
                                        periodicity = periodicity))) %>% 
    filter(rmse==min(rmse))
  
  out2 <- fit_model(data = data,cp = out1$cp,
                    model = model,
                    periodicity = periodicity,
                    return_fit = TRUE)
  
  out3 <- data %>% 
    mutate(before_cp = TRUE) %>%
    mutate(pred1=predict(out2$fit,newdata=.)) %>%
    mutate(before_cp = period>=out1$cp) %>% 
    mutate(pred2=predict(out2$fit,newdata=.)) %>% 
    select(period,n,pred1,pred2)
  
  return(list(cp=out1$cp,
              fit=out2$fit,
              pred=out3))
}

find_pred_bound_cp <- function(pred_data){
  pred_data %>% 
    arrange(period)  %>%
    mutate(above_cp=cumsum(n>pred1)==row_number()) %>%
    filter(above_cp) %>%
    summarise(cp=max(period)) %>% 
    .$cp
}

fit_cumsum_mods <- function(count_data, model, periodicity){
  cp_out <- count_data %>% 
    arrange(-period)
  
  t_series <- ts(cp_out$n, start = min(-1*cp_out$period), frequency = 1)
  
  cp_est <- suppressWarnings( cpts(cpt.mean(t_series,pen.value = 1, penalty = "None", test.stat = 'CUSUM')))
  cp <- cp_out$period[cp_est]
  
  fit <- fit_model(data = count_data,cp = cp,model = model,periodicity = periodicity,return_fit = TRUE)
  
  pred_data <-  count_data %>% 
    mutate(before_cp = TRUE) %>% 
    mutate(pred1 = predict(fit$fit, newdata = .))
  
  return(list(cp = cp,
              pred = pred_data))
}

get_rankings <- function(cp_data){
  tmp1 <- cp_data %>% 
    mutate(tmp=abs(cp_consistency)) %>% 
    distinct(tmp) %>% 
    arrange(tmp) %>% 
    mutate(consistency_rank=row_number()) %>% 
    inner_join(cp_data %>% 
                 mutate(tmp=abs(cp_consistency))) %>% 
    select(label,consistency_rank)
  
  tmp2 <- cp_data %>% 
    distinct(mse) %>% 
    arrange(mse) %>% 
    mutate(mse_rank=row_number()) %>% 
    inner_join(cp_data) %>% 
    select(label,mse_rank)
  
  tmp3 <- cp_data %>% 
    distinct(mse7) %>% 
    arrange(mse7) %>% 
    mutate(mse7_rank=row_number()) %>% 
    inner_join(cp_data) %>% 
    select(label,mse7_rank)
  
  tmp4 <- cp_data %>% 
    distinct(mse14) %>% 
    arrange(mse14) %>% 
    mutate(mse14_rank=row_number()) %>% 
    inner_join(cp_data) %>% 
    select(label,mse14_rank)
  
  cp_data %>% 
    select(label, cp, pb_cp = pred_bound_cp) %>% 
    left_join(tmp1, by = "label") %>% 
    left_join(tmp2, by = "label") %>% 
    left_join(tmp3, by = "label") %>% 
    left_join(tmp4, by = "label") %>% 
    mutate(avg_rank = (consistency_rank+mse_rank+mse7_rank+mse14_rank)/4) %>% 
    arrange(avg_rank) %>% 
    select(label,avg_rank,cp, pb_cp,consistency_rank:mse14_rank)
}


# count all visits
count_data_all <- visit_counts %>% 
  filter(is.na(dx_ver)) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))

# load(paste0(in_path,"/all_dx_visits.RData"))

#Change this to read from wherever the data is stored for you
count_data_ssd <- all_dx_visits %>%
  inner_join(ssd_codes, by = "dx") %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  left_join(tibble(days_since_index=min(visit_counts$days_since_index):-1),., by = "days_since_index") %>% # in case there are 0 days
  mutate(n = replace_na(n,0L)) %>% 
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))


# piecewise model results
ssd_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_ssd,
                                cp_range = 3:300,
                                model = .x,
                                periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))

# cusum model results
ssd_res2 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0("CUSUM ",label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~fit_cumsum_mods(count_data_ssd,model = .x, periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))

# aggregate results output
comb <- bind_rows(ssd_res1 %>% mutate(method = "piecewise"), 
                  ssd_res2 %>% mutate(method = "cusum"))

unique_cps <- distinct(comb, cp)

comb <- comb %>% 
  rename(true_cp = cp) %>% 
  mutate(cp = map(label,~tibble(cp=unique_cps$cp))) %>% 
  unnest(cp) %>% 
  filter(cp>=true_cp) %>% 
  mutate(mse = map2(pred,cp,
                    ~filter(.x, period>=.y) %>% # here its >=
                      summarise(mse = mean((n-pred1)^2)))) %>% 
  unnest(mse) %>% 
  select(-label) %>% 
  mutate(label = str_glue("Model = {model}, periodicity = {periodicity}, cp = {true_cp}, method = {method}")) 
  
mse_table <- comb %>% 
  arrange(desc(cp)) %>% 
  select(model, periodicity, true_cp, pred_bound_cp, method, cp, mse ) %>% 
  pivot_wider(id_cols = c("model", "periodicity",  "method","true_cp", "pred_bound_cp"), names_from = cp,  values_from = mse) %>% 
  rename(cp = true_cp) %>% 
  arrange(cp)

write_csv(mse_table, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", params$proj, "/tables/mse_table.csv"))


comb <- bind_rows(ssd_res1 %>% mutate(method = "piecewise"), 
                  ssd_res2 %>% mutate(method = "cusum"))

unique_cps <- distinct(comb, cp)

comb <- comb %>% 
  rename(true_cp = cp) %>% 
  mutate(cp = map(label,~tibble(cp=unique_cps$cp))) %>% 
  unnest(cp) %>% 
  filter(cp>=true_cp) %>% 
  mutate(mse = map2(pred,cp,
                    ~filter(.x, period>.y) %>% # here its > 
                      summarise(mse = mean((n-pred1)^2)))) %>% 
  unnest(mse) %>% 
  select(-label) %>% 
  mutate(label = str_glue("Model = {model}, periodicity = {periodicity}, cp = {true_cp}, method = {method}")) 

mse_table <- comb %>% 
  arrange(desc(cp)) %>% 
  select(model, periodicity, true_cp, pred_bound_cp, method, cp, mse ) %>% 
  pivot_wider(id_cols = c("model", "periodicity",  "method","true_cp", "pred_bound_cp"), names_from = cp,  values_from = mse) %>% 
  rename(cp = true_cp) %>% 
  arrange(cp)

write_csv(mse_table, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", params$proj, "/tables/mse_table_not_including_cp.csv"))