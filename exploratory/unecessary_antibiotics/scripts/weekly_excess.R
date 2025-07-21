

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

fit_model <- function(data,cp,model="lm"){
  
  tmp <- data %>% 
    filter(week>cp)
  
  if (model=="lm"){
    
    fit <- lm(n~week, data = tmp)
    
  } else if (model=="quad"){
    
    fit <- lm(n~poly(week,2), data = tmp)
    
  } else if (model=="cubic") {
    
    fit <- lm(n~poly(week,3), data = tmp)
    
  } else if (model=="exp"){
    
    fit <- lm(n~log(week), data = tmp)
  }
  
  pred <- data %>% 
    mutate(before_cp = TRUE) %>%
    mutate(pred=predict(fit,newdata=.))

    return(pred)
  }


###################
#### Load Data ####
###################

#### Load ABX visits -----------------------------------------------------------

# Note: this needs to be split off to better abx data location

abx_list <- read_xlsx("~/OneDrive - University of Iowa/delay_dx_projects/unecessary_antibiotics/data/antibiotics_for_delay_projectsPhil.xlsx") %>% 
  rename(name = `Antibiotic Name`) %>% 
  gather(key = disease, value = include, -name)

load("/Volumes/AML/truven_mind_projects/antibiotic_risk_categories/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == cond_name) %>% 
  filter(include == "Y") %>% 
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

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


#### Assemble Final Data -------------------------------------------------------

final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  distinct(patient_id,name,days_since_index) %>% 
    mutate(period = -days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) 

final_abx_counts <- final_abx_visits %>% 
  count(name,week)

# remove incomplete weeks
week_exclude <- final_abx_visits %>% 
  distinct(week,period) %>% 
  count(week) %>% 
  filter(n<7) %>% 
  select(week)

final_abx_counts <- final_abx_counts %>% 
  anti_join(week_exclude)

final_abx_visits <- final_abx_visits %>% 
  anti_join(week_exclude)

##############################
#### Fit Aggregate Models ####
##############################


## Fit models ------------------------------------------------------------------

abx_name_weekly_fits <- final_abx_counts %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(model = map(name,~c("lm","quad","cubic"))) %>% 
  unnest(model) %>% 
  ungroup() %>% 
  mutate(fits = map2(data,model,~fit_model(data = .x,
                                           cp = weekly_cp,
                                           model = .y))) %>% 
  select(name,model,fits) %>% 
  unnest(fits)

## Plots -----------------------------------------------------------------------

p1 <- final_abx_counts  %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw() +
  facet_wrap(~name, scale = "free_y") +
  ylab("Number of Prescriptions") +
  xlab("Weeks Before Diagnosis")


p2 <- abx_name_weekly_fits %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  geom_line(aes(y = pred, color = model), size = 1) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
  theme_bw() +
  facet_wrap(~name, scale = "free_y") +
  theme(legend.position = "bottom") +
  ylab("Number of Prescriptions") +
  xlab("Weeks Before Diagnosis")

p1
p2

p3 <- abx_name_weekly_fits %>% 
  filter(week <= weekly_cp) %>% 
  mutate(excess = n-pred) %>% 
  ggplot(aes(week,excess,color = model)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~name, scale = "free_y") +
  theme(legend.position = "bottom") +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Number of Excess Antibiotics") +
  xlab("Weeks Before Diagnosis")


#################################
#### Bootstrap Weekly trends ####
#################################



# simulation function
run_trial_weekly <- function(){
  tmp_index <- index_dates %>% 
    sample_frac(size = 1, replace = TRUE) %>% 
    mutate(boot_id = row_number())
  
  tmp_abx_counts <- final_abx_visits %>% 
    inner_join(tmp_index,by = join_by(patient_id),
               relationship = "many-to-many") %>% 
    distinct(boot_id,name,days_since_index,week) %>% 
    count(name,week)
  
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
                                             model = .y))) %>% 
    select(name,model,fits) %>% 
    unnest(fits) %>% 
    select(name,model,week,n,pred)
  
  # total abx fits
  total_fits <- tibble(model = c("lm","quad","cubic")) %>% 
    mutate(fits = map(model,~fit_model(data = tmp_abx_counts_total,
                                       cp = weekly_cp,
                                       model = .))) %>% 
    unnest(fits) %>% 
    select(model,week,n,pred) %>% 
    mutate(name = "total")
  
  tmp_fits <- bind_rows(named_fits,
                        total_fits)
  
  return(tmp_fits)
}

# run_trial_weekly() %>%
#   ggplot(aes(week,n)) +
#   geom_point() +
#   geom_line(aes(y=pred,color = model)) +
#   scale_x_reverse() +
#   facet_wrap(~name, scales = "free_y")

boot_res <- tibble(trial = 1:1000) %>% 
  mutate(res = map(trial,~run_trial_weekly())) %>% 
  unnest(res)

## Summarise Results Individual Antibiotics ------------------------------------

# start with individual and summed totals

# Note the total in this section reflects the individual results summed together

# mean fits and 95% upper bound
agg_boot_res <- boot_res %>% 
  group_by(model,week,name) %>% 
  summarise(mean_fit = mean(pred),
            hi_fit = quantile(pred,probs = 0.975),
            lo_fit = quantile(pred,probs = 0.025)) %>% 
  ungroup()

# individual abx
tmp1 <- agg_boot_res %>% 
  filter(name != "total") %>% 
  inner_join(final_abx_counts,by = join_by(week, name))

# individual sum totaled across models
tmp2 <- agg_boot_res %>% 
  filter(name != "total") %>% 
  group_by(model,week) %>%
  summarise(mean_fit = sum(mean_fit),
            hi_fit = sum(hi_fit),
            lo_fit = sum(lo_fit)) %>% 
  left_join(final_abx_counts %>% 
              group_by(week) %>% 
              summarise(n = sum(n))) %>% 
  mutate(name = "total (individual sum)") %>% 
  ungroup()

# modelled total
tmp3 <- agg_boot_res %>% 
  filter(name == "total") %>% 
  left_join(final_abx_counts %>% 
              group_by(week) %>% 
              summarise(n = sum(n)))

agg_boot_res <- bind_rows(tmp1,tmp2,tmp3)

agg_boot_res <- agg_boot_res %>% 
  mutate(excess = n-mean_fit,
         excess_high = n-lo_fit,
         excess_low = n-hi_fit) %>% 
  mutate_at(vars(excess:excess_low),~ifelse(week>weekly_cp,NA,.)) 

tmp1 <- agg_boot_res %>% 
  filter(week<=weekly_cp) %>% 
  group_by(model,name) %>% 
  summarise(excess = sum(excess),
            excess_low = sum(excess_low),
            excess_high = sum(excess_high),
            n = sum(n))

excess_res1 <- tmp1 %>% 
  mutate_at(vars(excess:excess_high),~round(.,2)) %>% 
  mutate(out = paste0(excess, " (", excess_low,", ",excess_high,")")) %>% 
  select(name,model,out) %>% 
  spread(key = model,value = out) %>% 
  select(abx = name, lm, quad, cubic) 

excess_res2 <- tmp1 %>% 
  mutate(frac = round(100*excess/n,2),
         frac_low = round(100*excess_low/n,2),
         frac_high = round(100*excess_high/n,2)) %>% 
  mutate_at(vars(frac:frac_high),~round(.,2)) %>% 
  mutate(out = paste0(frac, " (", frac_low,", ",frac_high,")")) %>% 
  select(name,model,out) %>% 
  spread(key = model,value = out) %>% 
  select(abx = name, lm, quad, cubic)

excess_res1
excess_res2


## plot model results by week --------------------------------------------------


plot_model_res <- function(model_name){
  
  plot_data <- filter(agg_boot_res,model == model_name)
  
  p1 <- plot_data %>% 
    mutate(delay_window = week<=weekly_cp) %>% 
    ggplot(aes(week,n)) +
    geom_point(aes(color = delay_window)) +
    scale_colour_manual(values = c("black","red")) +
    scale_x_reverse() +
    geom_line(aes(week, mean_fit), color = "blue") +
    geom_line(aes(week, hi_fit), color = "blue", linetype = 2) +
    geom_line(aes(week, lo_fit), color = "blue", linetype = 2) +
    facet_wrap(~name, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
    xlab("Weeks Before Diagnosis") +
    ylab("Number of Antibiotic Prescriptions")
  
  
  p2 <- plot_data %>% 
    mutate(delay_window = week<=weekly_cp) %>% 
    ggplot(aes(week,n)) +
    geom_point(aes(color = delay_window)) +
    scale_colour_manual(values = c("black","red")) +
    scale_x_reverse() +
    geom_ribbon(aes(ymin = lo_fit,ymax = hi_fit, alpha = 0.1)) +
    geom_line(aes(week, mean_fit), color = "blue") +
    facet_wrap(~name, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_vline(aes(xintercept = weekly_cp), linetype = 2) +
    xlab("Weeks Before Diagnosis") +
    ylab("Number of Antibiotic Prescriptions")
  
  p3 <- plot_data %>% 
    mutate(excess = n-mean_fit,
           excess_high = n-lo_fit,
           excess_low = n-hi_fit) %>% 
    filter(week<=weekly_cp) %>% 
    ggplot(aes(week,excess)) +
    geom_point() +
    geom_pointrange(aes(ymin = excess_low,ymax = excess_high)) +
    facet_wrap(~name) +
    scale_x_reverse() +
    theme_bw() +
    facet_wrap(~name, scales = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    xlab("Weeks Before Diagnosis") +
    ylab("Estimated Number of Excess Antibiotics")
  
  return(list(p1=p1,p2=p2,p3=p3))
}

lm_res_plots <- plot_model_res("lm")
quad_res_plots <- plot_model_res("quad")
cubic_res_plots <- plot_model_res("cubic")

lm_res_plots$p1
cubic_res_plots$p3

p4 <- agg_boot_res %>%  
  mutate(week = ifelse(model == "lm", week + 0.25,
                       ifelse(model == "cubic", week - 0.25, week))) %>% 
  filter(week<=weekly_cp) %>% 
  ggplot(aes(week,excess, color = model)) +
  geom_point() +
  geom_pointrange(aes(ymin = excess_low,
                      ymax = excess_high)) +
  scale_x_reverse() +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  ylab("Number of Excess Antibiotics") +
  xlab("Weeks Before Diagnosis")

p5 <- agg_boot_res %>% 
  mutate(frac = round(100*excess/n,2),
         frac_low = round(100*excess_low/n,2),
         frac_high = round(100*excess_high/n,2)) %>% 
  mutate(week = ifelse(model == "lm", week + 0.25,
                       ifelse(model == "cubic", week - 0.25, week))) %>% 
  filter(week<=weekly_cp) %>% 
  ggplot(aes(week,frac, color = model)) +
  geom_point() +
  geom_pointrange(aes(ymin = frac_low,
                      ymax = frac_high)) +
  scale_x_reverse() +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  ylab("% of Antibiotics Recieved in Excess") +
  xlab("Weeks Before Diagnosis")


################################
#### Simulate Excess Visits ####
################################

weekly_vis <- final_abx_visits %>% 
  select(patient_id,name,week) %>% 
  group_by(name,week) %>% 
  nest() %>% 
  ungroup()


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
  filter(week<=weekly_cp) %>% 
  filter(name!="total") %>% 
  select(trial,name,week,pred) %>% 
  inner_join(select(final_abx_counts,name,week,n),join_by(name, week)) %>% 
  mutate(excess = round(n-pred,0)) %>% 
  mutate(excess = ifelse(excess<0, 0, round(excess))) %>% 
  mutate(excess = ifelse(excess>n,n,excess)) %>% 
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

tmp_excess_lm <- tmp_excess


#### Simulate Excess Visits - Quadratic ----------------------------------------

# pull out lm results
tmp <- boot_res %>% 
  filter(model == "quad") %>% 
  filter(week<=weekly_cp) %>% 
  filter(name!="total") %>% 
  select(trial,name,week,pred) %>% 
  inner_join(select(final_abx_counts,name,week,n),join_by(name, week)) %>% 
  mutate(excess = round(n-pred,0)) %>% 
  mutate(excess = ifelse(excess<0, 0, round(excess))) %>% 
  mutate(excess = ifelse(excess>n,n,excess)) %>% 
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

quad_stats <- list(n_excess_patients = n_excess_patients,
                 frac_excess_patients = frac_excess_patients,
                 n_excess_patients_by_type = n_excess_patients_by_type,
                 frac_excess_patients_by_type = frac_excess_patients_by_type,
                 n_excess_abx_by_patient_count = n_excess_abx_by_patient_count,
                 mean_excess_per_patient_all = mean_excess_per_patient_all,
                 mean_excess_per_excess_patient = mean_excess_per_excess_patient)

tmp_excess_quad <- tmp_excess

#### Simulate Excess Visits - Cubic ----------------------------------------

# pull out lm results
tmp <- boot_res %>% 
  filter(model == "cubic") %>% 
  filter(week<=weekly_cp) %>% 
  filter(name!="total") %>% 
  select(trial,name,week,pred) %>% 
  inner_join(select(final_abx_counts,name,week,n),join_by(name, week)) %>% 
  mutate(excess = round(n-pred,0)) %>% 
  mutate(excess = ifelse(excess<0, 0, round(excess))) %>% 
  mutate(excess = ifelse(excess>n,n,excess)) %>% 
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

cubic_stats <- list(n_excess_patients = n_excess_patients,
                   frac_excess_patients = frac_excess_patients,
                   n_excess_patients_by_type = n_excess_patients_by_type,
                   frac_excess_patients_by_type = frac_excess_patients_by_type,
                   n_excess_abx_by_patient_count = n_excess_abx_by_patient_count,
                   mean_excess_per_patient_all = mean_excess_per_patient_all,
                   mean_excess_per_excess_patient = mean_excess_per_excess_patient)

tmp_excess_cubic <- tmp_excess


#### Save Results For Report ####

save(p1,p2,p3,
     lm_res_plots,quad_res_plots,cubic_res_plots,
     p4,p5,
     excess_res1,excess_res2,
     lm_stats,quad_stats,cubic_stats,
     file = "~/Data/projects/excess_abx/cocci/data/report_data.RData")

tibble(metric = c("Number of Patients with Excess Antibiotic",
                  "% of Patients with Excess Antibiotic",
                  "Excess Antibiotics per patient (all patients)",
                  "Excess Antibiotics per patient (patients w/ excess)"),
       
       lm = c(paste0(lm_stats$n_excess_patients$n_excess_patients," (",
                     lm_stats$n_excess_patients$n_excess_patients_lo,", ",
                     lm_stats$n_excess_patients$n_excess_patients_hi,")"),
              
              paste0(round(lm_stats$frac_excess_patients$n_excess_patients,2)," (",
                     round(lm_stats$frac_excess_patients$n_excess_patients_lo,2),", ",
                     round(lm_stats$frac_excess_patients$n_excess_patients_hi,2),")"),
              
              paste0(round(lm_stats$mean_excess_per_patient_all$mean_excess_per_patient,3)," (",
                     round(lm_stats$mean_excess_per_patient_all$mean_excess_per_patient_lo,3),", ",
                     round(lm_stats$mean_excess_per_patient_all$mean_excess_per_patient,3),")"),
              
              paste0(round(lm_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3)," (",
                     round(lm_stats$mean_excess_per_excess_patient$mean_excess_per_patient_lo,3),", ",
                     round(lm_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3),")")
       ),
       
       quad = c(paste0(quad_stats$n_excess_patients$n_excess_patients," (",
                       quad_stats$n_excess_patients$n_excess_patients_lo,", ",
                       quad_stats$n_excess_patients$n_excess_patients_hi,")"),
                
                paste0(round(quad_stats$frac_excess_patients$n_excess_patients,2)," (",
                       round(quad_stats$frac_excess_patients$n_excess_patients_lo,2),", ",
                       round(quad_stats$frac_excess_patients$n_excess_patients_hi,2),")"),
                
                paste0(round(quad_stats$mean_excess_per_patient_all$mean_excess_per_patient,3)," (",
                       round(quad_stats$mean_excess_per_patient_all$mean_excess_per_patient_lo,3),", ",
                       round(quad_stats$mean_excess_per_patient_all$mean_excess_per_patient,3),")"),
                
                paste0(round(quad_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3)," (",
                       round(quad_stats$mean_excess_per_excess_patient$mean_excess_per_patient_lo,3),", ",
                       round(quad_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3),")")
       ),
       
       cubic = c(paste0(cubic_stats$n_excess_patients$n_excess_patients," (",
                        cubic_stats$n_excess_patients$n_excess_patients_lo,", ",
                        cubic_stats$n_excess_patients$n_excess_patients_hi,")"),
                 
                 paste0(round(cubic_stats$frac_excess_patients$n_excess_patients,2)," (",
                        round(cubic_stats$frac_excess_patients$n_excess_patients_lo,2),", ",
                        round(cubic_stats$frac_excess_patients$n_excess_patients_hi,2),")"),
                 
                 paste0(round(cubic_stats$mean_excess_per_patient_all$mean_excess_per_patient,3)," (",
                        round(cubic_stats$mean_excess_per_patient_all$mean_excess_per_patient_lo,3),", ",
                        round(cubic_stats$mean_excess_per_patient_all$mean_excess_per_patient,3),")"),
                 
                 paste0(round(cubic_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3)," (",
                        round(cubic_stats$mean_excess_per_excess_patient$mean_excess_per_patient_lo,3),", ",
                        round(cubic_stats$mean_excess_per_excess_patient$mean_excess_per_patient,3),")")
                 
       ))


tmp_lm <- lm_stats$frac_excess_patients_by_type %>% 
  mutate_at(vars(n_excess_patients:n_excess_patients_hi),~round(.,2)) %>% 
  mutate(lm = paste0(n_excess_patients, " (",n_excess_patients_lo,", ",n_excess_patients_hi,")")) %>% 
  select(name, lm)


tmp_quad <- quad_stats$frac_excess_patients_by_type %>% 
  mutate_at(vars(n_excess_patients:n_excess_patients_hi),~round(.,2)) %>% 
  mutate(quad = paste0(n_excess_patients, " (",n_excess_patients_lo,", ",n_excess_patients_hi,")")) %>% 
  select(name, quad)

tmp_cubic <- cubic_stats$frac_excess_patients_by_type %>% 
  mutate_at(vars(n_excess_patients:n_excess_patients_hi),~round(.,2)) %>% 
  mutate(cubic = paste0(n_excess_patients, " (",n_excess_patients_lo,", ",n_excess_patients_hi,")")) %>% 
  select(name, cubic)

tmp_lm %>% 
  inner_join(tmp_quad, by = join_by(name)) %>% 
  inner_join(tmp_cubic, by = join_by(name))


tmp_lm <- lm_stats$n_excess_abx_by_patient_count %>% 
  mutate_at(vars(n_mean:n_hi),~round(.,2)) %>% 
  mutate(lm = paste0(n_mean, " (",n_lo,", ",n_hi,")")) %>% 
  select(n_excess, lm)

tmp_quad <- quad_stats$n_excess_abx_by_patient_count %>% 
  mutate_at(vars(n_mean:n_hi),~round(.,2)) %>% 
  mutate(quad = paste0(n_mean, " (",n_lo,", ",n_hi,")")) %>% 
  select(n_excess, quad)

tmp_cubic <- cubic_stats$n_excess_abx_by_patient_count %>% 
  mutate_at(vars(n_mean:n_hi),~round(.,2)) %>% 
  mutate(cubic = paste0(n_mean, " (",n_lo,", ",n_hi,")")) %>% 
  select(n_excess, cubic)

tmp_lm %>% 
  inner_join(tmp_quad, by = join_by(n_excess)) %>% 
  inner_join(tmp_cubic, by = join_by(n_excess))


### Add Days Supplied ####

tmp <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  distinct(patient_id,name,days_since_index,daysupp) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(week = 1+((period-1) %/% 7)) %>% 
  select(patient_id,name,week,daysupp)

tmp_excess$excess[[3]] %>% 
  inner_join(tmp) %>% 
  group_by(name) %>% 
  summarise(daysupp = sum(daysupp))



