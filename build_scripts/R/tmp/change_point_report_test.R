
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

#### peicewise models ----------------------------------------------------------
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


all_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_all,
                                cp_range = 3:300,
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


### Cusum models ---------------------------------------------------------------

ssd_res_out1 <- mutate(ssd_res1,
       tmp = ifelse(periodicity," w/ periodicity",""),
       label = paste0("Piecewise ", model, tmp)) %>% 
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
  select(label,cp,pred_bound_cp,mse,cp_n_miss,pb_cp_n_miss)

ssd_res_out1


ssd_res1$pred[[1]] %>% 
  filter(period<=58) %>% 
  mutate(., n_miss = ifelse(n>pred1,n-pred1,0)) %>% 
  summarise(sum(n_miss))


ssd_res1$pred[[1]] %>% 
  mutate(n_miss = ifelse(n>pred1,n-pred1,0)) 

library(trend)
install.packages("trend")

library(changepoint)
install.packages("changepoint")


cp_out <- count_data_ssd %>% 
  arrange(-period)

t_series <- ts(cp_out$n, start = min(-1*cp_out$period), frequency = 1)

cp_est <- cpts(cpt.mean(t_series,pen.value = 1, penalty = "None", test.stat = 'CUSUM'))
cp <- cp_out$period[cp_est]
cp

fit <- fit_model(data = count_data_ssd,cp = cp,model = "lm",periodicity = FALSE,return_fit = TRUE)

count_data_ssd %>% 
  mutate(before_cp = TRUE) %>% 
  mutate(pred1 = predict(fit$fit, newdata = .)) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_point() +
  geom_line(aes(y = pred1), color = "red") +
  geom_vline(aes(xintercept = -cp))

fit_cumsum_mods <- function(count_data, model, periodicity){
  cp_out <- count_data %>% 
    arrange(-period)
  
  t_series <- ts(cp_out$n, start = min(-1*cp_out$period), frequency = 1)
  
  cp_est <- suppressWarnings( cpts(cpt.mean(t_series,pen.value = 1, penalty = "None", test.stat = 'CUSUM')))
  cp <- cp_out$period[cp_est]
  
  fit <- fit_model(data = count_data_ssd,cp = cp,model = model,periodicity = periodicity,return_fit = TRUE)
  
  pred_data <-  count_data %>% 
    mutate(before_cp = TRUE) %>% 
    mutate(pred1 = predict(fit$fit, newdata = .))
  
  return(list(cp = cp,
              pred = pred_data))
}

fit_cumsum_mods(count_data_ssd,model = "lm",periodicity = TRUE)

tibble(model = c("lm","quad","cubic","exp")) %>% 
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




tmp3 <- ssd_res2 %>%
  filter(periodicity)

tmp4 <- ssd_res2 %>%
  filter(!periodicity)

tmp3 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp3, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp3, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()

tmp4 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp4, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp4, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()

tmp2 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred2), color = "blue") +
  geom_line(aes(y = pred1), color ="red") +
  geom_vline(data = tmp2, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp2, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()

#### SSD Visits ####

ssd_codes <- bind_rows(codeBuildr::load_ssd_codes(cond_name) %>%
                         filter(type == "icd9") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=9),
                       codeBuildr::load_ssd_codes(cond_name) %>%
                         filter(type == "icd10") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=10))


# load(paste0(in_path,"/all_dx_visits.RData"))

tibble(days_since_index=min(visit_counts$days_since_index):-1)

#Change this to read from wherever the data is stored for you
count_data_ssd <- all_dx_visits %>%
  inner_join(ssd_codes, by = "dx") %>%
  distinct(enrolid,days_since_index) %>%
  count(days_since_index) %>%
  left_join(tibble(days_since_index=min(visit_counts$days_since_index):-1),., by = "days_since_index") %>% 
  mutate(n = replace_na(n,0L)) %>% 
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))


tmp <- tibble(model = c("lm","quad","cubic","exp")) %>% 
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


tmp1 <- tmp %>%
  filter(periodicity)

tmp2 <- tmp %>%
  filter(!periodicity)
