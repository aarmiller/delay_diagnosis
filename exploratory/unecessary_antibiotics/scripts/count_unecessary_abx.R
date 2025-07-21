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

boot_res <- tibble(trial = 1:1000) %>% 
  mutate(res = map(trial,~run_trial_weekly()))

boot_res <- boot_res %>% 
  unnest(res)

#### Summarise Results ####


# mean fits and 95% upper bound
agg_boot_res <- boot_res %>% 
  group_by(model,period,name) %>% 
  summarise(mean_fit = mean(pred),
            hi_fit = quantile(pred,probs = 0.975),
            lo_fit = quantile(pred,probs = 0.025)) %>% 
  ungroup()

# count of excess antibiotics
excess_res1 <- agg_boot_res %>% 
  left_join(abx_counts_name_weekly) %>% 
  filter(period<=weekly_cp) %>% 
  mutate(excess = n-mean_fit,
            excess_high = n-lo_fit,
            excess_low = n-hi_fit) %>% 
  group_by(model,name) %>% 
  summarise(excess = sum(excess),
            excess_low = sum(excess_low),
            excess_high = sum(excess_high)) %>% 
  mutate_at(vars(excess:excess_high),~round(.,2)) %>% 
  mutate(out = paste0(excess, " (", excess_low,", ",excess_high,")")) %>% 
  select(name,model,out) %>% 
  spread(key = model,value = out) %>% 
  select(abx = name, lm, quad, cubic) %>% 
  filter(abx!="total")

# percent of antibiotics before that were in excess
tmp1 <- agg_boot_res %>% 
  left_join(abx_counts_name_weekly) %>% 
  filter(period<=weekly_cp) %>% 
  mutate(excess = n-mean_fit,
         excess_high = n-lo_fit,
         excess_low = n-hi_fit) %>% 
  group_by(model,name) %>% 
  summarise(excess = sum(excess),
            excess_low = sum(excess_low),
            excess_high = sum(excess_high)) %>% 
  ungroup()

excess_res2 <- abx_counts_name_weekly %>% 
  filter(period<=weekly_cp) %>% 
  group_by(name) %>% 
  summarise(total_abx = sum(n)) %>% 
  inner_join(tmp1) %>% 
  gather(key = key,value = value, -name, -total_abx, -model) %>% 
  mutate(frac = round(100*value/total_abx,2)) %>% 
  select(-value) %>% 
  spread(key = key, value =frac) %>% 
  mutate(out = paste0(excess, " (", excess_low,", ",excess_high,")")) %>% 
  select(name,model,out) %>% 
  spread(key = model,value = out) %>% 
  select(abx = name, lm, quad, cubic)

excess_res1
excess_res2


## linear model results --------------------------------------------------------
tmp1 <- bind_rows(abx_counts_name_weekly,
                  mutate(abx_counts_total_weekly,name = "total")) %>% 
  inner_join(filter(agg_boot_res,model == "lm"))

tmp2 <- select(filter(boot_res,model == "lm"),trial,name,period,pred) %>% 
  inner_join(tmp1) 

bind_rows(abx_counts_name_weekly,
          mutate(abx_counts_total_weekly,name = "total")) %>% 
  inner_join(filter(agg_boot_res,model == "lm")) %>% 
  mutate(delay_window = period<=weekly_cp) %>% 
  ggplot(aes(period,n)) +
  geom_point(aes(color = delay_window)) +
  scale_colour_manual(values = c("black","red")) +
  scale_x_reverse() +
  geom_line(aes(period, mean_fit), color = "blue") +
  geom_line(aes(period, hi_fit), color = "blue", linetype = 2) +
  geom_line(aes(period, lo_fit), color = "blue", linetype = 2) +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
  xlab("Weeks Before Diagnosis") +
  ylab("Number of Antibiotic Prescriptions")
  

bind_rows(abx_counts_name_weekly) %>% 
  inner_join(filter(agg_boot_res,model == "lm")) %>% 
  mutate(delay_window = period<=weekly_cp) %>% 
  ggplot(aes(period,n)) +
  geom_point(aes(color = delay_window)) +
  scale_colour_manual(values = c("black","red")) +
  scale_x_reverse() +
  geom_ribbon(aes(ymin = lo_fit,ymax = hi_fit, alpha = 0.1)) +
  geom_line(aes(period, mean_fit), color = "blue") +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
  xlab("Weeks Before Diagnosis") +
  ylab("Number of Antibiotic Prescriptions")

bind_rows(abx_counts_name_weekly) %>% 
  inner_join(filter(agg_boot_res,model == "lm")) %>% 
  mutate(excess = n-mean_fit,
         excess_high = n-lo_fit,
         excess_low = n-hi_fit) %>% 
  filter(period<=weekly_cp) %>% 
  ggplot(aes(period,excess)) +
  geom_point() +
  geom_pointrange(aes(ymin = excess_low,ymax = excess_high)) +
  facet_wrap(~name) +
  scale_x_reverse() +
  theme_bw() +
  facet_wrap(~name, scales = "free_y") +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  xlab("Weeks Before Diagnosis") +
  ylab("Estimated Number of Excess Antibiotics")

################################
#### Simulate Excess Visits ####
################################

weekly_vis <- final_abx_visits %>% 
  select(patient_id,name,week) %>% 
  group_by(name,week) %>% 
  nest() %>% 
  ungroup()

abx_counts_name_weekly %>% 
  filter(name == "metronidazole",
         week == 3)

final_abx_visits %>% 
  filter(name == "metronidazole",
         week == 3)

draw_excess <- function(excess_counts){
  weekly_vis %>% 
    inner_join(excess_counts,by = join_by(name, week)) %>% 
    mutate(excess = map2(data,excess,~sample_n(.x, size =.y, replace = FALSE))) %>% 
    select(name,week,excess) %>% 
    unnest(excess)
}

compute_excess_stats <- function(excess_draw_data){
  tmp1 <- excess_draw_data %>% 
    summarise(n_excess_abx_patients = n_distinct(patient_id))
  
  tmp2 <- excess_draw_data %>% 
    group_by(name) %>% 
    summarise(n_excess_patients = n_distinct(patient_id))
  
  tmp3 <- index_dates %>% 
    distinct(patient_id) %>% 
    left_join(excess_draw_data %>% 
                count(patient_id,name = "n_excess"),
              by = join_by(patient_id)) %>% 
    mutate(n_excess = replace_na(n_excess,0L)) %>% 
    count(n_excess)
  
  tmp4 <- index_dates %>% 
    distinct(patient_id) %>% 
    left_join(excess_draw_data %>% 
                count(patient_id,name = "n_excess"),
              by = join_by(patient_id)) %>% 
    mutate(n_excess = replace_na(n_excess,0L)) %>% 
    summarise(mean_excess = mean(n_excess))
  
  tmp5 <- excess_draw_data %>% 
    count(patient_id,name = "n_excess") %>% 
    summarise(mean_excess = mean(n_excess))
  
  
  out <- list(n_excess_abx_patients = tmp1,
              n_excess_patients_by_type = tmp2,
              n_excess_abx_by_patient = tmp3,
              mean_excess_per_patient_all = tmp4,
              mean_excess_per_excess_patient = tmp5)
  
  return(out)
}

#### Simulate Excess Visits - Linear -------------------------------------------

# pull out lm results
tmp <- boot_res %>% 
  filter(model == "lm") %>% 
  filter(period<=weekly_cp) %>% 
  filter(name!="total") %>% 
  select(trial,name,week = period,pred) %>% 
  inner_join(select(abx_counts_name_weekly,name,week,n),join_by(name, week)) %>% 
  mutate(excess = round(n-pred,0)) %>% 
  mutate(excess = ifelse(excess<0, 0, round(excess))) %>% 
  select(trial,week,name,excess) %>% 
  group_by(trial) %>% 
  nest()

# compute excess stats
tmp_excess <- tmp %>% 
  mutate(excess = map(data,draw_excess)) %>% 
  mutate(excess_stats = map(excess,compute_excess_stats))



## aggregate final counts

# number of patients with excess abx
n_excess_patients <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_patients)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(n_excess_patients = mean(n_excess_abx_patients),
            n_excess_patients_lo = quantile(n_excess_abx_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_abx_patients, probs = 0.975))

# fraction of patients with excess abx
frac_excess_patients <- n_excess_patients %>% 
  mutate_all(~100*(./nrow(index_dates)))

# number of patients with excess abx by type
n_excess_patients_by_type <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_patients_by_type)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(n_excess_patients_lo = quantile(n_excess_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_patients, probs = 0.975),
            n_excess_patients = mean(n_excess_patients)) %>% 
  select(name,n_excess_patients,everything())

# fraction of patients with excess abx by type
frac_excess_patients_by_type <- n_excess_patients_by_type %>% 
  mutate_at(vars(n_excess_patients:n_excess_patients_hi),~100*(./nrow(index_dates)))

# count of number of excess antibiotics by patient
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_by_patient)) %>% 
  select(trial,res) %>% 
  unnest(res)

tmp_span <-  select(tmp_excess,trial) %>% 
  mutate(n_excess = map(trial,~0:max(tmp$n_excess))) %>% 
  unnest(n_excess)

n_excess_abx_by_patient_count <- tmp_span %>% 
  left_join(tmp) %>% 
  mutate(n = replace_na(n,0L)) %>% 
  group_by(n_excess) %>% 
  summarise(n_mean = mean(n),
            n_lo = quantile(n,probs = 0.025),
            n_hi = quantile(n,probs = 0.975))

tmp_excess

# mean number of excess antibiotics by patient across all patients
mean_excess_per_patient_all <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_patient_all)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(mean_excess_per_patient = mean(mean_excess),
            mean_excess_per_patient_lo = quantile(mean_excess, probs = 0.025),
            mean_excess_per_patient_hi = quantile(mean_excess, probs = 0.975))

# mean number of excess antibiotics by patient across patients with excess
mean_excess_per_excess_patient <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_excess_patient)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(mean_excess_per_patient = mean(mean_excess),
            mean_excess_per_patient_lo = quantile(mean_excess, probs = 0.025),
            mean_excess_per_patient_hi = quantile(mean_excess, probs = 0.975))

## Aggregate stats

lm_stats <- list(n_excess_patients = n_excess_patients,
                 frac_excess_patients = frac_excess_patients,
                 n_excess_patients_by_type = n_excess_patients_by_type,
                 frac_excess_patients_by_type = frac_excess_patients_by_type,
                 n_excess_abx_by_patient_count = n_excess_abx_by_patient_count,
                 mean_excess_per_patient_all = mean_excess_per_patient_all,
                 mean_excess_per_excess_patient = mean_excess_per_excess_patient)


#### Simulate Excess Visits - Quadratic ----------------------------------------

# pull out lm results
tmp <- boot_res %>% 
  filter(model == "quad") %>% 
  filter(period<=weekly_cp) %>% 
  filter(name!="total") %>% 
  select(trial,name,week = period,pred) %>% 
  inner_join(select(abx_counts_name_weekly,name,week,n),join_by(name, week)) %>%
  mutate(excess = ifelse(pred<0,n,round(n-pred))) %>% 
  mutate(excess = ifelse(excess<0, 0, excess)) %>% 
  select(trial,week,name,excess) %>% 
  group_by(trial) %>% 
  nest()

# compute excess stats
tmp_excess <- tmp %>% 
  mutate(excess = map(data,draw_excess)) %>% 
  mutate(excess_stats = map(excess,compute_excess_stats))

tmp %>% 
  unnest(data) %>% 
  ungroup() %>% 
  inner_join(tmp2) %>% 
  filter(excess > avail_draw_count)

tmp2 <- weekly_vis %>% 
  mutate(avail_draw_count = map_int(data,~nrow(.))) %>% 
  select(name,week,avail_draw_count)


select(abx_counts_name_weekly,name,week,n) %>% 
  filter(week == 3)

weekly_vis %>% 
  filter(name == "metronidazole") %>% 
  filter(week == 3)




## aggregate final counts

# number of patients with excess abx
n_excess_patients <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_patients)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(n_excess_patients = mean(n_excess_abx_patients),
            n_excess_patients_lo = quantile(n_excess_abx_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_abx_patients, probs = 0.975))

# fraction of patients with excess abx
frac_excess_patients <- n_excess_patients %>% 
  mutate_all(~100*(./nrow(index_dates)))

# number of patients with excess abx by type
n_excess_patients_by_type <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_patients_by_type)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  group_by(name) %>% 
  summarise(n_excess_patients_lo = quantile(n_excess_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_patients, probs = 0.975),
            n_excess_patients = mean(n_excess_patients)) %>% 
  select(name,n_excess_patients,everything())

# fraction of patients with excess abx by type
frac_excess_patients_by_type <- n_excess_patients_by_type %>% 
  mutate_at(vars(n_excess_patients:n_excess_patients_hi),~100*(./nrow(index_dates)))

# count of number of excess antibiotics by patient
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_by_patient)) %>% 
  select(trial,res) %>% 
  unnest(res)

tmp_span <-  select(tmp_excess,trial) %>% 
  mutate(n_excess = map(trial,~0:max(tmp$n_excess))) %>% 
  unnest(n_excess)

n_excess_abx_by_patient_count <- tmp_span %>% 
  left_join(tmp) %>% 
  mutate(n = replace_na(n,0L)) %>% 
  group_by(n_excess) %>% 
  summarise(n_mean = mean(n),
            n_lo = quantile(n,probs = 0.025),
            n_hi = quantile(n,probs = 0.975))

tmp_excess

# mean number of excess antibiotics by patient across all patients
mean_excess_per_patient_all <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_patient_all)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(mean_excess_per_patient = mean(mean_excess),
            mean_excess_per_patient_lo = quantile(mean_excess, probs = 0.025),
            mean_excess_per_patient_hi = quantile(mean_excess, probs = 0.975))

# mean number of excess antibiotics by patient across patients with excess
mean_excess_per_excess_patient <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_excess_patient)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(mean_excess_per_patient = mean(mean_excess),
            mean_excess_per_patient_lo = quantile(mean_excess, probs = 0.025),
            mean_excess_per_patient_hi = quantile(mean_excess, probs = 0.975))

## Aggregate stats

lm_stats <- list(n_excess_patients = n_excess_patients,
                 frac_excess_patients = frac_excess_patients,
                 n_excess_patients_by_type = n_excess_patients_by_type,
                 frac_excess_patients_by_type = frac_excess_patients_by_type,
                 n_excess_abx_by_patient_count = n_excess_abx_by_patient_count,
                 mean_excess_per_patient_all = mean_excess_per_patient_all,
                 mean_excess_per_excess_patient = mean_excess_per_excess_patient)


