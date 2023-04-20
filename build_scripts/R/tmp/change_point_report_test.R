
rm(list = ls())
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "tb"

load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

in_path <- paste0(delay_params$path,"delay_results")
out_path <- paste0(delay_params$path,"change_point_results/")

if (!dir.exists(out_path)){
  dir.create(out_path)
}


### Set bounds on cp analysis
cp_lower <- ifelse(is.na(delay_params$cp_lower),3,delay_params$cp_lower)
cp_upper <- ifelse(is.na(delay_params$cp_upper),round(.80*delay_params$upper_bound,0),delay_params$cp_upper)

######################
##### Get Counts #####
######################

# load("/Volumes/AML/small_dbs/tb/truven/enroll_restrict_365/delay_results/all_dx_visits.RData")

# count all visits
count_data_all <- visit_counts %>% 
  filter(is.na(dx_ver)) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))

# load SSD codes
ssd_codes <- bind_rows(codeBuildr::load_ssd_codes(cond_name) %>%
                         filter(type == "icd9") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=9),
                       codeBuildr::load_ssd_codes(cond_name) %>%
                         filter(type == "icd10") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=10))


# load(paste0(in_path,"/all_dx_visits.RData"))

#Change this to read from wherever the data is stored for you
count_data_ssd <- all_dx_visits %>%
  inner_join(ssd_codes, by = "dx") %>%
  distinct(enrolid,days_since_index) %>%
  count(days_since_index) %>%
  left_join(tibble(days_since_index=min(visit_counts$days_since_index):-1),., by = "days_since_index") %>% # in case there are 0 days
  mutate(n = replace_na(n,0L)) %>% 
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))

###################
#### Functions ####
###################


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

# tmp <- find_cp(count_data,
#                cp_range = 3:300,
#                model = "exp",
#                periodicity = FALSE)

##################################################
#### Get Change-Point Results for all visits #####
##################################################

#### ssd visits ----------------------------------------------------------------
ssd_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_ssd,
                                cp_range = cp_lower:cp_upper,
                                model = .x,
                                periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))

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

#### all visits ----------------------------------------------------------------
all_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_all,
                                cp_range = cp_lower:cp_upper,
                                model = .x,
                                periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))

all_res2 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0("CUSUM ",label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~fit_cumsum_mods(count_data_all,model = .x, periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))

###########################
#### Aggregate Results ####
###########################

ssd_res_out <- bind_rows(mutate(ssd_res1,
                                 tmp = ifelse(periodicity," w/ periodicity",""),
                                 label = paste0("Piecewise ", model, tmp)),
                          mutate(ssd_res2,
                                 tmp = ifelse(periodicity," w/ periodicity",""),
                                 label = paste0("CUSUM ", model, tmp))) %>% 
  select(label,cp,pred_bound_cp,pred) %>% 
  mutate(pred = map(pred,~mutate(., n_miss = ifelse(n>pred1,n-pred1,0)))) %>% 
  mutate(cp_n_miss = map2(pred,cp,~filter(.x, period<=.y) %>% summarise(cp_n_miss = sum(n_miss)))) %>% 
  unnest(cp_n_miss) %>% 
  mutate(pb_cp_n_miss = map2(pred,pred_bound_cp,~filter(.x, period<=.y) %>% summarise(pb_cp_n_miss = sum(n_miss)))) %>% 
  unnest(pb_cp_n_miss) %>% 
  mutate(mse = map2(pred,cp,
                    ~filter(.x, period>.y) %>% 
                      summarise(mse = mean((n-pred1)^2)))) %>%
  unnest(mse) %>% 
  mutate(mse7 = map2(pred,cp,
                    ~filter(.x, between(period,.y,.y+6)) %>% 
                      summarise(mse7 = mean((n-pred1)^2)))) %>% 
  unnest(mse7) %>% 
  mutate(mse14 = map2(pred,cp,
                     ~filter(.x, between(period,.y,.y+13)) %>% 
                       summarise(mse14 = mean((n-pred1)^2)))) %>% 
  unnest(mse14) %>% 
  select(label,cp,pred_bound_cp,mse,cp_n_miss,pb_cp_n_miss,mse7,mse14) %>% 
  mutate_at(vars(mse:mse14),~round(.,2)) %>% 
  mutate(cp_consistency = 100*(cp-pred_bound_cp)/cp)

all_res_out <- bind_rows(mutate(all_res1,
                                tmp = ifelse(periodicity," w/ periodicity",""),
                                label = paste0("Piecewise ", model, tmp)),
                         mutate(all_res2,
                                tmp = ifelse(periodicity," w/ periodicity",""),
                                label = paste0("CUSUM ", model, tmp))) %>% 
  select(label,cp,pred_bound_cp,pred) %>% 
  mutate(pred = map(pred,~mutate(., n_miss = ifelse(n>pred1,n-pred1,0)))) %>% 
  mutate(cp_n_miss = map2(pred,cp,~filter(.x, period<=.y) %>% summarise(cp_n_miss = sum(n_miss)))) %>% 
  unnest(cp_n_miss) %>% 
  mutate(pb_cp_n_miss = map2(pred,pred_bound_cp,~filter(.x, period<=.y) %>% summarise(pb_cp_n_miss = sum(n_miss)))) %>% 
  unnest(pb_cp_n_miss) %>% 
  mutate(mse = map2(pred,cp,
                    ~filter(.x, period>.y) %>% 
                      summarise(mse = mean((n-pred1)^2)))) %>% 
  unnest(mse) %>% 
  mutate(mse7 = map2(pred,cp,
                     ~filter(.x, between(period,.y,.y+6)) %>% 
                       summarise(mse7 = mean((n-pred1)^2)))) %>% 
  unnest(mse7) %>% 
  mutate(mse14 = map2(pred,cp,
                      ~filter(.x, between(period,.y,.y+13)) %>% 
                        summarise(mse14 = mean((n-pred1)^2)))) %>% 
  unnest(mse14) %>% 
  select(label,cp,pred_bound_cp,mse,cp_n_miss,pb_cp_n_miss,mse7,mse14) %>% 
  mutate_at(vars(mse:mse14),~round(.,2)) %>% 
  mutate(cp_consistency = 100*(cp-pred_bound_cp)/cp)

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

ssd_res_out %>% 
  select(label,cp,pred_bound_cp,cp_n_miss,pb_cp_n_miss)

get_rankings(ssd_res_out)
get_rankings(all_res_out)

ssd_res_out %>% 
  distinct(mse7) %>% 
  arrange(mse7) %>% 
  mutate(mse7_rank=row_number()) %>% 
  inner_join(ssd_res_out) %>% 
  select(label,mse7_rank)

all_res_out %>% 
  mutate(cp_consistency = 100*(cp-pred_bound_cp)/cp)

opt_mods <- get_rankings(ssd_res_out) %>% 
  filter(avg_rank == min(avg_rank))

opt_mods_text <- ifelse(nrow(opt_mods)>1,
                            paste0("The optimal models were ",paste0(opt_mods$label,collapse = " or ")),
                            paste0("The optimal model was ",opt_mods$label))

paste0("The optimal models were ",paste0(opt_mods$label,collapse = " or "))

opt_mods %>% 
  select(label) %>% 
  inner_join(ssd_res_out) %>% 
  select(label,`Change Point`=cp,`Pred. Bound CP` = pred_bound_cp,`CP # Miss` = cp_n_miss,`PB CP # Miss` = pb_cp_n_miss)

tmp <- bind_rows(mutate(ssd_res1,
                        tmp = ifelse(periodicity," w/ periodicity",""),
                        label = paste0("Piecewise ", model, tmp)),
                 mutate(ssd_res2,
                        tmp = ifelse(periodicity," w/ periodicity",""),
                        label = paste0("CUSUM ", model, tmp))) %>% 
  inner_join(select(opt_mods,label)) %>% 
  select(label,cp,pred_bound_cp,pred)

tmp %>% 
  unnest(pred) %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()
