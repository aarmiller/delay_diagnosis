rm(list = ls())

library(tidyverse)
library(readxl)


################
#### Params ####
################


cond_name <- "cocci"

cp <- 91

weekly_cp <- round(cp/7,0)

load("/Volumes/AML/params/delay_any_params.RData")

upper_bound <- delay_any_params[[cond_name]]$upper_bound

###################
#### Functions ####
###################

fit_model <- function(data,cp,model="lm",periodicity = FALSE, return_fit=FALSE, return_pred = FALSE){
  
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
  
  if (!return_fit & !return_pred){
    return(rmse)
  } else {
    
    out <- list(rmse = rmse)
    
    if (return_fit){
      out[["fit"]] <- fit
    }
    
    if (return_pred){
      
      pred <- data %>% 
        mutate(before_cp = TRUE) %>%
        mutate(pred1=predict(fit,newdata=.)) %>%
        mutate(before_cp = period>=cp) %>% 
        mutate(pred2=predict(fit,newdata=.)) %>% 
        select(period,n,pred1,pred2)
      
      out[["pred"]] <- pred
    }
    
    return(out)
    
  }
}

###################
#### Load Data ####
###################

#### Load ABX visits -----------------------------------------------------------

# Note: this needs to be split off to better abx data location

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx")

load("/Volumes/AML/truven_mind_projects/antibiotic_risk_categories/antibiotics_groupings_new.RData")

ndc_codes <- antibiotic_ndc_groups_new$ndcnum


#### Load Disease Data ---------------------------------------------------------

db <- src_sqlite(paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

abx_visits <- db %>% 
  tbl("all_rx_visits") %>% 
  filter(ndcnum %in% ndc_codes) %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

index_dates <- index_dates %>% 
  filter(time_before_index>=upper_bound) %>% 
  distinct(patient_id,index_date)


####################################
#### Compute Overall ABX Counts ####
####################################

## Compute Counts --------------------------------------------------------------

abx_counts_total <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(select(filter(abx_list,cocci=="Y"),name = `Antibiotic Name`)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  distinct(patient_id,ndcnum,name,days_since_index) %>% 
  count(days_since_index) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

abx_counts_name <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(select(filter(abx_list,cocci=="Y"),name = `Antibiotic Name`)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  distinct(patient_id,ndcnum,name,days_since_index) %>% 
  count(name,days_since_index) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

## Fit models ------------------------------------------------------------------

abx_total_daily_fits <- tibble(model = c("lm","quad","cubic")) %>% 
  mutate(fits = map(model,~fit_model(data = abx_counts_total,
                                     cp = cp,
                                     model = .,
                                     periodicity = TRUE,
                                     return_pred = TRUE)$pred)) %>% 
  unnest(fits) 

abx_name_daily_fits <- abx_counts_name %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(model = map(name,~c("lm","quad","cubic"))) %>% 
  unnest(model) %>% 
  ungroup() %>% 
  mutate(fits = map2(data,model,~fit_model(data = .x,
                                           cp = cp,
                                           model = .y,
                                           periodicity = TRUE,
                                           return_pred = TRUE)$pred)) %>% 
  select(name,model,fits) %>% 
  unnest(fits)


#### Plots ---------------------------------------------------------------------

# Aggregate Counts
p1 <- abx_counts_total %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw()

# Aggregate Counts with fits
p2 <- abx_total_daily_fits %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = pred1, color = model), size = 1) +
  theme_bw() +
  geom_vline(aes(xintercept = cp), linetype = 2)

# Named Counts
p3 <- abx_counts_name %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw() +
  facet_wrap(~name, scales = "free_y") 

# Named Counts with fits
p4 <- abx_name_daily_fits %>% 
  ggplot(aes(period,n)) +
  geom_point(alpha = 0.25) +
  geom_line(aes(y = pred1, color = model),size = 1) +
  scale_x_reverse() +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  geom_vline(aes(xintercept = cp))

p1 
p2
p3
p4

###################################
#### Compute Weekly ABX Counts ####
###################################

## Get counts ------------------------------------------------------------------

## Exclude incomplete weeks
week_exclude <- abx_counts_total %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)
  

abx_counts_total_weekly <- abx_counts_total %>% 
  anti_join(week_exclude, by = "week") %>% 
  group_by(week) %>% 
  summarise(n = sum(n)) %>% 
  mutate(period = week)

abx_counts_name_weekly <- abx_counts_name %>% 
  anti_join(week_exclude, by = "week") %>% 
  group_by(name,week) %>% 
  summarise(n = sum(n)) %>% 
  ungroup() %>% 
  mutate(period = week)

## Fit models ------------------------------------------------------------------

abx_total_weekly_fits <- tibble(model = c("lm","quad","cubic")) %>% 
  mutate(fits = map(model,~fit_model(data = abx_counts_total_weekly,
                                     cp = weekly_cp,
                                     model = .,
                                     periodicity = FALSE,
                                     return_pred = TRUE)$pred)) %>% 
  unnest(fits) 

abx_name_weekly_fits <- abx_counts_name_weekly %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(model = map(name,~c("lm","quad","cubic"))) %>% 
  unnest(model) %>% 
  ungroup() %>% 
  mutate(fits = map2(data,model,~fit_model(data = .x,
                                           cp = weekly_cp,
                                           model = .y,
                                           periodicity = FALSE,
                                           return_pred = TRUE)$pred)) %>% 
  select(name,model,fits) %>% 
  unnest(fits)



## Plots -----------------------------------------------------------------------

p5 <- abx_counts_total_weekly %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw()

p6 <- abx_total_weekly_fits %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred1, color = model)) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
  theme_bw() +
  xlab("Weeks Before Diagnosis") +
  ylab("Number of relevant antibiotic prescribed")

p6_1 <- abx_total_weekly_fits %>% 
  mutate(excess_abx = n-pred1) %>% 
  filter(period<=weekly_cp) %>% 
  ggplot(aes(period,excess_abx)) +
  geom_line(aes(color = model)) +
  scale_x_reverse() +
  theme_bw() +
  xlab("Weeks Before Diagnosis") +
  ylab("Excess Antibiotics")

p7 <- abx_name_weekly_fits %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw() +
  facet_wrap(~name, scale = "free_y")

p8 <- abx_name_weekly_fits %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred1, color = model), size = 1) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
  theme_bw() +
  facet_wrap(~name, scale = "free_y") +
  theme(legend.position = "bottom")

p6
p6_1

tmp <- abx_name_weekly_fits %>% 
  mutate(excess = n-pred1) %>% 
  mutate(excess = ifelse(excess<0,0,excess)) %>% 
  filter(period<=weekly_cp) %>% 
  group_by(period,model) %>% 
  mutate(excess_frac = 100*excess/sum(excess)) %>% 
  mutate(abx_frac = 100*n/sum(n)) %>% 
  ungroup()

bind_rows(select(tmp,name,model,period,frac = excess_frac) %>% 
            mutate(group = "Excess ABX Fraction"),
          select(tmp,name,model,period,frac = abx_frac) %>% 
            mutate(group = "Fraction of ABX Prescribed")) %>% 
  filter(model == "lm") %>% 
  ggplot(aes(period,frac,color = group)) +
  geom_point() +
  facet_wrap(~name, scale = "free_y") +
  theme_bw() +
  scale_x_reverse() +
  theme(legend.position = "bottom") +
  ylab("Fraction of Antibiotics Drawn")


#################################
#### Bootstrap Weekly trends ####
#################################

final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(select(filter(abx_list,cocci=="Y"),name = `Antibiotic Name`)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  distinct(patient_id,name,days_since_index) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  mutate(week = 1+((period-1) %/% 7)) %>% 
  anti_join(week_exclude)

rm(tmp_index,tmp_abx_counts,tmp_abx_counts_total,named_fits)

# simulation function
run_trial_weekly <- function(){
  tmp_index <- index_dates %>% 
    sample_frac(size = 1, replace = TRUE) %>% 
    mutate(boot_id = row_number())
  
  tmp_abx_counts <- final_abx_visits %>% 
    inner_join(tmp_index,by = join_by(patient_id),
               relationship = "many-to-many") %>% 
    distinct(boot_id,name,days_since_index,week) %>% 
    count(name,week) %>% 
    mutate(period = week)
  
  tmp_abx_counts_total <- tmp_abx_counts %>% 
    group_by(week) %>% 
    summarise(n = sum(n)) %>% 
    mutate(period = week)
  
  # named abx fits
  
  named_fits <- tmp_abx_counts %>% 
    group_by(name) %>% 
    nest() %>% 
    mutate(model = map(name,~c("lm","quad","cubic"))) %>% 
    unnest(model) %>% 
    ungroup() %>% 
    mutate(fits = map2(data,model,~fit_model(data = .x,
                                             cp = weekly_cp,
                                             model = .y,
                                             periodicity = FALSE,
                                             return_pred = TRUE)$pred)) %>% 
    select(name,model,fits) %>% 
    unnest(fits) %>% 
    select(model,period,n,pred = pred1,name)
  
  # total abx fits
  total_fits <- tibble(model = c("lm","quad","cubic")) %>% 
    mutate(fits = map(model,~fit_model(data = tmp_abx_counts_total,
                                       cp = weekly_cp,
                                       model = .,
                                       periodicity = FALSE,
                                       return_pred = TRUE)$pred)) %>% 
    unnest(fits) %>% 
    select(model,period,n,pred = pred1) %>% 
    mutate(name = "total")

  tmp_fits <- bind_rows(total_fits,
                        named_fits)
  
  return(tmp_fits)
}

# run_trial_weekly() %>% 
#   ggplot(aes(period,n)) +
#   geom_point() +
#   geom_line(aes(y=pred,color = model)) +
#   scale_x_reverse() +
#   facet_wrap(~name, scales = "free_y")

boot_res <- tibble(trial = 1:500) %>% 
  mutate(res = map(trial,~run_trial_weekly()))

boot_res <- boot_res %>% 
  unnest(res)

#### Summarise Results ####


# mean fits and 95% upper bound
agg_boot_res <- boot_res %>% 
  group_by(model,period,name) %>% 
  summarise(mean_fit = mean(pred),
            hi95_fit = quantile(pred,probs = 0.95)) %>% 
  ungroup()

## linear model results --------------------------------------------------------
tmp1 <- bind_rows(abx_counts_name_weekly,
                  mutate(abx_counts_total_weekly,name = "total")) %>% 
  inner_join(filter(agg_boot_res,model == "lm"))

tmp2 <- select(filter(boot_res,model == "lm"),trial,name,period,pred) %>% 
  inner_join(tmp1) 

tmp2 %>% 
  mutate(excess = ifelse(n>hi95_fit & week <= weekly_cp,
                         TRUE,FALSE)) %>% 
  ggplot(aes(period,n)) +
  geom_point(aes(color = excess)) +
  scale_colour_manual(values = c("black","red")) +
  scale_x_reverse() +
  geom_line(aes(period, pred, group = trial), color = "blue", alpha = 0.1) +
  geom_line(aes(period, hi95_fit), color = "red") +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") 

bind_rows(abx_counts_name_weekly,
          mutate(abx_counts_total_weekly,name = "total")) %>% 
  inner_join(filter(agg_boot_res,model == "lm")) %>% 
  mutate(excess = ifelse(n>hi95_fit & week <= weekly_cp,
                         TRUE,FALSE)) %>% 
  ggplot(aes(period,n)) +
  geom_point(aes(color = excess)) +
  scale_colour_manual(values = c("black","red")) +
  scale_x_reverse() +
  geom_line(aes(period, hi95_fit), color = "red") +
  geom_line(aes(period, mean_fit), color = "blue") +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2)
  

tmp1 %>% 
  mutate(excess_hi = n-hi95_fit,
         excess_mean = n-mean_fit) %>% 
  filter(period<=weekly_cp) %>% 
  group_by(name) %>% 
  summarise(excess_hi = sum(excess_hi),
            excess_mean = sum(excess_mean))

# proportion of each antibiotic that is in excess
# note this shows why we can simply draw at random
tmp1 %>% 
  mutate(excess_hi = n-hi95_fit,
         excess_mean = n-mean_fit) %>% 
  filter(period<=weekly_cp) %>% 
  mutate(excess_frac = excess_mean/n) %>% 
  ggplot(aes(period,excess_frac)) +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  scale_x_reverse() +
  theme_bw() +
  xlab("Weeks Before Diagnosis") +
  ylab("Fraction of prescriptions in excess")

# Fraction of antibiotics prescribed by class
tmp1 %>% 
  filter(name != "total") %>% 
  group_by(week) %>% 
  mutate(total = sum(n)) %>% 
  mutate(abx_frac = 100*n/total) %>% 
  ggplot(aes(period,abx_frac)) +
  geom_line() +
  scale_x_reverse() +
  facet_wrap(~name) +
  ylab("Fraction of Antibiotics Prescribed") +
  theme_bw() +
  xlab("Weeks Before Diagnosis")

## What we should be drawing

tmp1 %>% 
  mutate(excess_hi = n-hi95_fit,
         excess_mean = n-mean_fit) %>%
  filter(period<=weekly_cp)  %>% 
  select(name,week,period,n,excess_mean) %>% 
  filter(name != "total") %>% 
  group_by(week) %>% 
  mutate(draw_total = sum(excess_mean)) %>% 
  mutate(frac_draw = 100*excess_mean/draw_total) %>% 
  ggplot(aes(period,frac_draw)) +
  geom_point() +
  facet_wrap(~name) +
  scale_x_reverse() +
  theme_bw() +
  xlab("Weeks Before Diagnosis") +
  ylab("Fraction of prescriptions in excess")


bind_rows(tmp1 %>% 
            filter(name != "total") %>% 
            group_by(week) %>% 
            mutate(total = sum(n)) %>% 
            mutate(frac = 100*n/total) %>% 
            select(name,period,frac) %>% 
            ungroup() %>% 
            mutate(group = "abx_frac"),
          tmp1 %>% 
            mutate(excess_hi = n-hi95_fit,
                   excess_mean = n-mean_fit) %>%
            filter(period<=weekly_cp)  %>% 
            select(name,week,period,n,excess_mean,excess_hi) %>% 
            filter(name != "total") %>% 
            group_by(week) %>% 
            # mutate(draw_total = sum(excess_mean)) %>% 
            # mutate(frac_draw = 100*excess_mean/draw_total) %>% 
            mutate(draw_total = sum(excess_hi)) %>% 
            mutate(frac_draw = 100*excess_hi/draw_total) %>% 
            select(name,period,frac = frac_draw) %>% 
            mutate(group = "draw_frac")) %>% 
  filter(period<=weekly_cp) %>% 
  ggplot(aes(period,frac,color =group)) +
  geom_point() +
  scale_x_reverse() +
  facet_wrap(~name, scale = "free_y") +
  theme_bw()



############ Exploratory ###########


tmp_excess <- tmp1 %>% 
  mutate(excess_hi = n-hi95_fit,
         excess_mean = n-mean_fit) %>% 
  filter(period<=weekly_cp) %>% 
  select(name,week,excess = excess_mean) %>% 
  mutate(excess = ifelse(excess<0, 0, round(excess)))


sim_visits <- function(excess_counts){
  tmp <- final_abx_visits %>% 
    group_by(name,week) %>% 
    nest() %>% 
    ungroup() %>% 
    inner_join(excess_counts,by = join_by(name, week)) %>% 
    mutate(data = map2(data,excess,~sample_n(.x,size = .y,replace = FALSE))) %>% 
    select(name,week,data) %>% 
    unnest(data) 
  
  out1 <- tmp %>% 
    summarise(n_pat = n_distinct(patient_id),
              n_abx = n())
  
  out2 <- tmp %>%  
    group_by(name) %>% 
    summarise(n_pat = n_distinct(patient_id),
              n_abx = n())
  
  return(list(out1,
              out2))
}

sim_visits(tmp_excess)


sim_vis_res <- tibble(trial = 1:100) %>% 
  mutate(res = map(trial,~sim_visits(tmp_excess)))

sim_vis_res <- sim_vis_res %>% 
  mutate(out1 = map(res,~.[[1]])) %>% 
  mutate(out2 = map(res,~.[[2]])) %>% 
  select(trial,out1,out2)

sim_vis_res %>% 
  select(trial,out1) %>% 
  unnest(out1) %>% 
  summarise(n_pat_mean = mean(n_pat),
            n_pat_lo = quantile(n_pat,probs = 0.025),
            n_pat_hi = quantile(n_pat,probs = 0.975))



## Compute excess using expected from each trial --------

tmp2 %>% 
  filter(week<=weekly_cp) %>% 
  mutate(excess = n-pred) %>% 
  group_by(trial,name) %>% 
  summarise(excess = sum(excess)) %>% 
  group_by(name) %>% 
  summarise(mean_excess = mean(excess),
            lo_excess = quantile(excess,probs = 0.025),
            high_excess = quantile(excess,probs = 0.975))



