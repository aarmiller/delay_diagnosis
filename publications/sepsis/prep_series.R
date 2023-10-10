
library(tidyverse)

rm(list = ls())

load("/Volumes/AML/params/final_delay_params.RData")

cond_name <- "sarcoid_skin"

delay_params <- final_delay_params[[cond_name]]

ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>% 
  filter(type %in% c("icd9","icd10")) %>% 
  mutate(dx_ver = ifelse(type=="icd9",9L,10L)) %>% 
  select(dx = code,dx_ver)

load(paste0("/Volumes/Statepi_Diagnosis/projects/sarcoid/",cond_name,"/index_cases.RData"))
load(paste0("/Volumes/Statepi_Diagnosis/prelim_results/sarcoid/delay_results/all_dx_visits.RData"))

index_cases

patient_ids <- index_cases %>% distinct(patient_id)
all_dx_visits <- all_dx_visits %>% inner_join(patient_ids)

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

tmp <- all_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index)

visit_counts <- bind_rows(tmp,visit_counts %>%
                            filter(days_since_index<=0))

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


count_data_ssd %>% 
  ggplot(aes(days_since_index,n)) +
  geom_point()

bound_val <- 730-500

rm(cp1,tmp_data,fit_lm,fit_cube,fit_quad,fit_exp,pred_data_out,
   pred_data_tmp,find_pred_bound_cp,implied_cps,find_fit_measures,
   out1)

cp1 <- 153

fit_trends <- function(count_data,cp1){
  
  # filter to training data
  tmp_data <- count_data_ssd %>% 
    filter(period>=cp1)
  
  # fit models
  fit_lm <- lm(n~period+dow, data = tmp_data)
  fit_quad <- lm(n~poly(period,2)+dow, data = tmp_data)
  fit_cube <- lm(n~poly(period,3)+dow, data = tmp_data)
  fit_exp <- lm(n~log(period)+dow, data = tmp_data)
  
  # get fitted values
  pred_data_tmp <- count_data_ssd %>% 
    mutate(lm = predict(fit_lm,newdata =.),
           quad = predict(fit_quad,newdata =.),
           cube = predict(fit_cube,newdata =.),
           exp = predict(fit_exp,newdata =.)) 
  
  find_pred_bound_cp <- function(pred_data){
    tmp_lm <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>lm)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_quad <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>quad)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_cube <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>cube)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_exp <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>exp)==row_number()) %>% 
      filter(above_cp) 
    
    tibble(model = c("lm","quad","cube","exp"),
           implied_cp = c(ifelse(nrow(tmp_lm)==0,NA,max(tmp_lm$period)),
                          ifelse(nrow(tmp_quad)==0,NA,max(tmp_quad$period)),
                          ifelse(nrow(tmp_cube)==0,NA,max(tmp_cube$period)),
                          ifelse(nrow(tmp_exp)==0,NA,max(tmp_exp$period))))
  }
  
  implied_cps <- find_pred_bound_cp(pred_data_tmp)
  
  pred_data_out <- pred_data_tmp %>% 
    select(period,n,lm:exp) %>% 
    gather(key = "model", value = "pred", -period, -n) %>% 
    inner_join(implied_cps, by = "model") 
  
  
  find_fit_measures <- function(pred_data,bound_val){
    pred_data %>%
      group_by(model) %>%
      filter(period>=implied_cp) %>% 
      summarise(mse=mean((n-pred)^2))
  }
  
  out1 <- inner_join(implied_cps,find_fit_measures(pred_data_out,cp1), by = "model")
  
  return(list(cp_res=out1,
              pred = pred_data_out))
  
}

fit_trends(count_data_ssd,121)

tmp <- tibble(cp = 50:400) %>% 
  mutate(res = map(cp,~fit_trends(count_data_ssd,.)))

tmp <- tmp %>% 
  mutate(cp_res = map(res,~.$cp_res),
         pred = map(res,~.$pred)) %>% 
  select(-res)


find_cp <- function(data,cp_deviation,cp_fp){
  
  tmp <- data %>% 
    select(cp,cp_res) %>% 
    unnest(cp_res) %>% 
    mutate(cp_diff = abs(cp-implied_cp)) %>% 
    filter(implied_cp<=(cp-cp_deviation),
           cp_diff<=cp_fp)            # allow implied change-point to deviate to the left of the provided cp
    
    if (nrow(tmp)>0){ 
      return(group_by(tmp,model) %>% 
               filter(mse==min(mse)))
    } else {
        tmp
      }
}



tmp_res <- tibble(cp_deviation = 1:10) %>% 
  mutate(cp_fp=map(cp_deviation,~1:15)) %>% 
  unnest(cp_fp) %>% 
  mutate(res = map2(cp_deviation,cp_fp,~find_cp(tmp,.x,.y)))

tmp_res %>% 
  unnest(res) %>% 
  ggplot(aes(cp_deviation,cp_fp, fill = mse)) +
  geom_tile() +
  geom_text(aes(label = implied_cp)) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~model)


tmp_res %>% 
  unnest(res) %>% 
  filter(model=="pred_exp")



tmp %>% 
  select(cp,cp_res) %>% 
  unnest(cp_res) %>% 
  mutate(cp_diff = abs(cp-implied_cp)) %>% 
  filter(implied_cp<=(cp-1)) %>%            # allow implied change-point to deviate to the left of the provided cp
  filter(cp_diff<=1) %>%                   # allow the cp and implied cp to differ by given amount
  group_by(model) %>% 
  filter(mse==min(mse))
  arrange(mse)
  


tmp2 <- fit_trends(count_data_ssd,383) 

tmp2$pred %>% 
  select(period,n,pred_lm:pred_exp) %>% 
  gather(key = "model", value = "pred", -period, -n) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = pred, color = model))

