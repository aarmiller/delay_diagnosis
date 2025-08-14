

rm(list = ls())

library(tidyverse)
library(readxl)
library(smallDB)

################
#### Params ####
################

cond_name <- "histo"

weekly_cp <- 15

upper_bound <- 365

# load("/Volumes/AML/params/delay_any_params.RData")

# upper_bound <- delay_any_params[[cond_name]]$upper_bound


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

load("/Volumes/argon_home/projects/delay_diagnosis/excess_abx/data/antibiotics_groupings_new.RData")

abx_include <- abx_list %>% 
  filter(disease == cond_name) %>% 
  filter(include %in% c("Y","?")) %>% 
  filter(!(name %in% c("cefaclor","cefadroxil","cefditoren","cefixime",
                       "cefpodoxime"))) %>% 
  filter(!(name %in% c("amoxicillin","cefprozil","cefuroxime","cephalexin","clindamycin"))) %>%  # remove due to negative CI
  select(name)

ndc_codes <- antibiotic_ndc_groups_new %>% 
  inner_join(abx_include) %>% 
  .$ndcnum

#### Load Disease Data ---------------------------------------------------------

db <- DBI::dbConnect(RSQLite::SQLite(), paste0("~/Data/MarketScan/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

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

crosswalk <- db %>% 
  tbl("enrolid_crosswalk") %>% 
  collect()

histo_cw <- index_dates %>% 
  inner_join(crosswalk)

save(histo_cw, file = "/Volumes/AML/tmp_transfer/excess_abx/histo_index_cases.RData")

## Pull Medicaid RX enroll On Argon --------------------------------------------

load("/Shared/AML/tmp_transfer/excess_abx/histo_index_cases.RData")

db <- DBI::dbConnect(RSQLite::SQLite(),"/Shared/Statepi_Marketscan/databases/Truven/enrollment_dbs/medicaid_fixes.db")

tmp_ids <- histo_cw %>% 
  filter(medicaid==1) %>% 
  select(enrolid)

tmp_collect_enroll_med <- function(year){
  db %>% 
    tbl(paste0("medicaid_drugcovg_",year)) %>% 
    inner_join(tmp_ids, copy = TRUE, by = join_by(enrolid)) %>% 
    collect()
}

tmp_enroll <- tibble(year=14:21) %>% 
  mutate(data = map(year,tmp_collect_enroll_med))

med_rx_enroll <- tmp_enroll %>% 
  unnest(data) %>% 
  inner_join(select(histo_cw,patient_id,enrolid)) %>% 
  select(-year)

save(med_rx_enroll, file = "/Shared/AML/tmp_transfer/excess_abx/histo_med_rx_enroll.RData")

## Pull rx enroll locally for ccae/mdcr ----------------------------------------

tmp_collect_enroll <- function(source,year){
  db %>% 
    tbl(paste0("enrollment_detail_",source,"_",year)) %>% 
    filter(rx==1) %>% 
    select(patient_id,dtstart,dtend) %>% 
    collect() 
}

tmp_enroll_info <- collect_plan(db) %>% 
  filter(source!="medicaid") %>% 
  mutate(data = map2(source,year,tmp_collect_enroll))

load("/Volumes/AML/tmp_transfer/excess_abx/histo_med_rx_enroll.RData")


# assemble rx enrollment dates
rx_enroll <- bind_rows(tmp_enroll_info %>% 
                         unnest(data) %>% 
                         select(patient_id,dtstart,dtend),
                       select(med_rx_enroll,patient_id,dtstart,dtend)) %>% 
  arrange(patient_id,dtstart) %>% 
  mutate(gap = dtstart != lag(dtend)+1) %>% 
  mutate(gap = replace_na(gap,FALSE)) %>% 
  mutate(new_id = patient_id!=lag(patient_id)) %>% 
  mutate(new_id = replace_na(new_id,FALSE)) %>% 
  mutate(tmp =  ifelse(gap==TRUE | new_id == TRUE, 1L, 0L)) %>% 
  mutate(period = cumsum(tmp)) %>% 
  group_by(patient_id,period) %>% 
  summarise(dtstart = min(dtstart),
            dtend = max(dtend)) %>% 
  ungroup()

# Final index dates
index_dates <- index_dates %>% 
  inner_join(rx_enroll) %>% 
  filter(index_date<=dtend & index_date>=dtstart) %>% 
  mutate(rx_days_before = index_date-dtstart) %>% 
  filter(rx_days_before>=365) %>% 
  select(patient_id,index_date) 


#### Assemble Final Data -------------------------------------------------------

# distinct abx visit days (to use for analysis)
final_abx_visits <- abx_visits %>% 
  inner_join(select(antibiotic_ndc_groups_new,ndcnum,name)) %>% 
  inner_join(index_dates) %>% 
  mutate(days_since_index = date - index_date) %>% 
  filter(between(days_since_index,-upper_bound,-1)) %>% 
  group_by(patient_id,name,days_since_index) %>% 
  summarise(daysupp = sum(daysupp)) %>% 
  ungroup() %>% 
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

################################
#### Compute Baseline Stats ####
################################


demo_data <- bind_rows(db %>% 
                         tbl("all_enroll") %>% 
                         select(patient_id,dobyr,sex,enrmon) %>% 
                         collect(),
                       db %>% 
                         tbl("all_enroll_medicaid") %>% 
                         select(patient_id,dobyr,sex,enrmon) %>% 
                         collect())

demo_data <- demo_data %>% 
  group_by(patient_id) %>% 
  filter(enrmon == max(enrmon)) %>% 
  ungroup()

demo_data <- index_dates %>% 
  inner_join(demo_data) %>% 
  mutate(age = year(as_date(index_date))-dobyr) %>% 
  mutate(age_cat = cut(age,breaks = c(-1,17,34,49,64,120)))

tmp1 <- demo_data %>% 
  summarise(out = as.character(n())) %>% 
  mutate(name = "N",
         cat = "") %>% 
  select(name,cat,out)

tmp2 <- demo_data %>% 
  count(sex) %>% 
  mutate(frac = round(100*n/sum(n),2)) %>% 
  mutate(cat = ifelse(sex==1,"Male","Female"),
         name = "sex") %>% 
  mutate(out = paste0(n," (",frac,")")) %>% 
  select(name,cat,out)

tmp3 <- demo_data %>% 
  summarise(mean_age = round(mean(age),2),
            median_age = median(age)) %>% 
  mutate(name = "age",
         cat = "mean/median",
         out = paste0(mean_age," (",median_age,")")) %>% 
  select(name,cat,out)

tmp4 <- demo_data %>% 
  count(age_cat) %>% 
  mutate(frac = round(100*n/sum(n),2)) %>% 
  mutate(out = paste0(n," (",frac,")")) %>% 
  mutate(name = "age") %>% 
  select(name,cat=age_cat,out)


tmp5 <- db %>% 
  tbl("tm") %>% 
  select(patient_id,index_date=svcdate,mdcr:medicaid) %>% 
  collect() %>% 
  inner_join(index_dates) %>% 
  summarise(mdcr = sum(mdcr),
            ccae = sum(ccae),
            medicaid = sum(medicaid)) %>% 
  gather(key = cat, value = n) %>% 
  mutate(name = "source") %>% 
  mutate(frac = round(100*n / sum(n),2)) %>% 
  mutate(out = paste0(n," (",frac,")")) %>% 
  select(name,cat,out)


tmp_abx_count <- filter(final_abx_visits,
                        week<=15) %>% 
  count(patient_id, name = "abx_count")

tmp6 <- demo_data %>% 
  left_join(tmp_abx_count) %>% 
  mutate(abx_count = replace_na(abx_count,0)) %>% 
  summarise(abx = sum(abx_count>0),
            abx_frac = round(100*abx/n(),2)) %>% 
  mutate(name = "abx",
         cat = "1 - 15 weeks",
         out = paste0(abx," (", abx_frac,")")) %>% 
  select(name,cat,out)

tmp_abx_count <- filter(final_abx_visits,
                        week<=52) %>% 
  count(patient_id, name = "abx_count")

tmp7 <- demo_data %>% 
  left_join(tmp_abx_count) %>% 
  mutate(abx_count = replace_na(abx_count,0)) %>% 
  summarise(abx = sum(abx_count>0),
            abx_frac = round(100*abx/n(),2)) %>% 
  mutate(name = "abx",
         cat = "1 - 52 weeks",
         out = paste0(abx," (", abx_frac,")")) %>% 
  select(name,cat,out)

tmp_abx_count <- filter(final_abx_visits,
                        week<=52 & week>15) %>% 
  count(patient_id, name = "abx_count")

tmp8 <- demo_data %>% 
  left_join(tmp_abx_count) %>% 
  mutate(abx_count = replace_na(abx_count,0)) %>% 
  summarise(abx = sum(abx_count>0),
            abx_frac = round(100*abx/n(),2)) %>% 
  mutate(name = "abx",
         cat = "16 - 52 weeks",
         out = paste0(abx," (", abx_frac,")")) %>% 
  select(name,cat,out)


baseline_data <- bind_rows(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)

save(baseline_data, file = paste0("~/Data/projects/excess_abx/",cond_name,"/data/baseline_data.RData"))

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
lm_res_plots$p2
lm_res_plots$p3


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

p4
p5

################################
#### Simulate Excess Visits ####
################################

weekly_vis <- final_abx_visits %>% 
  select(patient_id,name,week,daysupp) %>% 
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
    summarise(n_excess_abx_patients = n_distinct(patient_id),
              n_excess_abx_daysupp = sum(daysupp))
  
  tmp2 <- excess_draw_data %>% 
    group_by(name) %>% 
    summarise(n_excess_patients = n_distinct(patient_id),
              n_excess_daysupp = sum(daysupp))
  
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

# draw_excess(tmp$data[[1]]) %>% 
#   compute_excess_stats()

# compute excess stats
tmp_excess <- tmp %>% 
  mutate(excess = map(data,draw_excess)) %>% 
  mutate(excess_stats = map(excess,compute_excess_stats))

## aggregate final counts

# number of patients with excess abx
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_patients)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>%
  mutate(excess_daysupp_pp = n_excess_abx_daysupp/n_excess_abx_patients,
         excess_daysupp_pp_all = n_excess_abx_daysupp/nrow(index_dates)) %>%  
  summarise(n_excess_patients_mean = mean(n_excess_abx_patients),
            n_excess_patients_median = median(n_excess_abx_patients),
            n_excess_patients_lo = quantile(n_excess_abx_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_abx_patients, probs = 0.975),
            
            tot_excess_daysupp_mean = mean(n_excess_abx_daysupp),
            tot_excess_daysupp_median = median(n_excess_abx_daysupp),
            tot_excess_daysupp_lo = quantile(n_excess_abx_daysupp, probs = 0.025),
            tot_excess_daysupp_hi = quantile(n_excess_abx_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975))

# number of patients with excess abx
excess_patients <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure, 
         mean = n_excess_patients_mean,
         median = n_excess_patients_median,
         lo = n_excess_patients_lo,
         hi = n_excess_patients_hi) 

# fraction of patients with excess abx
frac_excess_patients <- excess_patients %>% 
  mutate(measure = "Percent of Patients with Excess ABX") %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates)))

# number of excess days supplied
excess_daysupp <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure, 
         mean = tot_excess_daysupp_mean,
         median = tot_excess_daysupp_median,
         lo = tot_excess_daysupp_lo,
         hi = tot_excess_daysupp_hi)

# number of excess days supplied per person
excess_daysupp_pp <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

# number of excess days supplied per person
excess_daysupp_pp_all <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# compute excess numbers by antibiotic type
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_patients_by_type)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  mutate(excess_daysupp_pp = n_excess_daysupp/n_excess_patients,
         excess_daysupp_pp_all = n_excess_daysupp/nrow(index_dates)) %>% 
  group_by(name) %>% 
  summarise(excess_patients_mean = mean(n_excess_patients),
            excess_patients_median = median(n_excess_patients),
            excess_patients_lo = quantile(n_excess_patients, probs = 0.025),
            excess_patients_hi = quantile(n_excess_patients, probs = 0.975),
            
            excess_daysupp_mean = mean(n_excess_daysupp),
            excess_daysupp_median = median(n_excess_daysupp),
            excess_daysupp_lo = quantile(n_excess_daysupp, probs = 0.025),
            excess_daysupp_hi = quantile(n_excess_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975)) 

# number of excess patients by type
excess_patients_by_type <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure,
         name,
         mean = excess_patients_mean,
         median = excess_patients_median,
         lo = excess_patients_lo,
         hi = excess_patients_hi)

# fraction of patients with excess abx by type
frac_excess_patients_by_type <- excess_patients_by_type %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates))) %>% 
  mutate(measure = "Percent of Patients with Excess ABX")

# excess daysupplied by type
excess_daysupp_by_type <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure,
         name,
         mean = excess_daysupp_mean,
         median = excess_daysupp_median,
         lo = excess_daysupp_lo,
         hi = excess_daysupp_hi)

# excess days supplied per person by type
excess_daysupp_pp_by_type <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

excess_daysupp_pp_all_by_type <- tmp %>%
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# bind_rows(excess_patients_by_type,
#           frac_excess_patients_by_type,
#           excess_daysupp_by_type,
#           excess_daysupp_pp_by_type,
#           excess_daysupp_pp_all_by_type) %>% 
#   mutate_at(vars(mean,lo,hi),~round(.,2)) %>% 
#   mutate(out = paste0(mean, " (",lo,", ",hi,")")) %>% 
#   select(measure,name,out) %>% 
#   spread(key = measure, value = out) %>% 
#   select(name,
#          `Number of Patients with Excess ABX`,
#          `Percent of Patients with Excess ABX`,
#          `Total Excess Days Supplied`,
#          `Excess Days Supplied Per Patient (excess patients)`,
#          `Excess Days Supplied Per Patient (all patients)`)

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
  summarise(mean = mean(n),
            median = median(n),
            lo = quantile(n,probs = 0.025),
            hi = quantile(n,probs = 0.975))


# mean number of excess antibiotics by patient across all patients
mean_excess_per_patient_all <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_patient_all)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (all patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

# mean number of excess antibiotics by patient across patients with excess
mean_excess_per_excess_patient <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_excess_patient)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (excess patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

## Aggregate stats
n_excess_abx_by_patient_count

lm_stats <- list(main_excess_stats = bind_rows(excess_patients,
                                               frac_excess_patients,
                                               excess_daysupp,
                                               excess_daysupp_pp,
                                               excess_daysupp_pp_all,
                                               mean_excess_per_excess_patient,
                                               mean_excess_per_patient_all),
                 
                 main_excess_stats_by_type = bind_rows(excess_patients_by_type,
                                                       frac_excess_patients_by_type,
                                                       excess_daysupp_by_type,
                                                       excess_daysupp_pp_by_type,
                                                       excess_daysupp_pp_all_by_type),
                 
                 excess_distribution = n_excess_abx_by_patient_count)

tmp_excess_lm <- tmp_excess


#### Simulate Excess Visits - Quadratic -------------------------------------------

# pull out quadratic results
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
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_patients)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>%
  mutate(excess_daysupp_pp = n_excess_abx_daysupp/n_excess_abx_patients,
         excess_daysupp_pp_all = n_excess_abx_daysupp/nrow(index_dates)) %>%  
  summarise(n_excess_patients_mean = mean(n_excess_abx_patients),
            n_excess_patients_median = median(n_excess_abx_patients),
            n_excess_patients_lo = quantile(n_excess_abx_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_abx_patients, probs = 0.975),
            
            tot_excess_daysupp_mean = mean(n_excess_abx_daysupp),
            tot_excess_daysupp_median = median(n_excess_abx_daysupp),
            tot_excess_daysupp_lo = quantile(n_excess_abx_daysupp, probs = 0.025),
            tot_excess_daysupp_hi = quantile(n_excess_abx_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975))

# number of patients with excess abx
excess_patients <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure, 
         mean = n_excess_patients_mean,
         median = n_excess_patients_median,
         lo = n_excess_patients_lo,
         hi = n_excess_patients_hi) 

# fraction of patients with excess abx
frac_excess_patients <- excess_patients %>% 
  mutate(measure = "Percent of Patients with Excess ABX") %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates)))

# number of excess days supplied
excess_daysupp <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure, 
         mean = tot_excess_daysupp_mean,
         median = tot_excess_daysupp_median,
         lo = tot_excess_daysupp_lo,
         hi = tot_excess_daysupp_hi)

# number of excess days supplied per person
excess_daysupp_pp <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

# number of excess days supplied per person
excess_daysupp_pp_all <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# compute excess numbers by antibiotic type
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_patients_by_type)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  mutate(excess_daysupp_pp = n_excess_daysupp/n_excess_patients,
         excess_daysupp_pp_all = n_excess_daysupp/nrow(index_dates)) %>% 
  group_by(name) %>% 
  summarise(excess_patients_mean = mean(n_excess_patients),
            excess_patients_median = median(n_excess_patients),
            excess_patients_lo = quantile(n_excess_patients, probs = 0.025),
            excess_patients_hi = quantile(n_excess_patients, probs = 0.975),
            
            excess_daysupp_mean = mean(n_excess_daysupp),
            excess_daysupp_median = median(n_excess_daysupp),
            excess_daysupp_lo = quantile(n_excess_daysupp, probs = 0.025),
            excess_daysupp_hi = quantile(n_excess_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975)) 

# number of excess patients by type
excess_patients_by_type <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure,
         name,
         mean = excess_patients_mean,
         median = excess_patients_median,
         lo = excess_patients_lo,
         hi = excess_patients_hi)

# fraction of patients with excess abx by type
frac_excess_patients_by_type <- excess_patients_by_type %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates))) %>% 
  mutate(measure = "Percent of Patients with Excess ABX")

# excess daysupplied by type
excess_daysupp_by_type <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure,
         name,
         mean = excess_daysupp_mean,
         median = excess_daysupp_median,
         lo = excess_daysupp_lo,
         hi = excess_daysupp_hi)

# excess days supplied per person by type
excess_daysupp_pp_by_type <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

excess_daysupp_pp_all_by_type <- tmp %>%
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# bind_rows(excess_patients_by_type,
#           frac_excess_patients_by_type,
#           excess_daysupp_by_type,
#           excess_daysupp_pp_by_type,
#           excess_daysupp_pp_all_by_type) %>% 
#   mutate_at(vars(mean,lo,hi),~round(.,2)) %>% 
#   mutate(out = paste0(mean, " (",lo,", ",hi,")")) %>% 
#   select(measure,name,out) %>% 
#   spread(key = measure, value = out) %>% 
#   select(name,
#          `Number of Patients with Excess ABX`,
#          `Percent of Patients with Excess ABX`,
#          `Total Excess Days Supplied`,
#          `Excess Days Supplied Per Patient (excess patients)`,
#          `Excess Days Supplied Per Patient (all patients)`)

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
  summarise(mean = mean(n),
            median = median(n),
            lo = quantile(n,probs = 0.025),
            hi = quantile(n,probs = 0.975))


# mean number of excess antibiotics by patient across all patients
mean_excess_per_patient_all <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_patient_all)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (all patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

# mean number of excess antibiotics by patient across patients with excess
mean_excess_per_excess_patient <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_excess_patient)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (excess patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

## Aggregate stats
quad_stats <- list(main_excess_stats = bind_rows(excess_patients,
                                               frac_excess_patients,
                                               excess_daysupp,
                                               excess_daysupp_pp,
                                               excess_daysupp_pp_all,
                                               mean_excess_per_excess_patient,
                                               mean_excess_per_patient_all),
                 
                 main_excess_stats_by_type = bind_rows(excess_patients_by_type,
                                                       frac_excess_patients_by_type,
                                                       excess_daysupp_by_type,
                                                       excess_daysupp_pp_by_type,
                                                       excess_daysupp_pp_all_by_type),
                 
                 excess_distribution = n_excess_abx_by_patient_count)

tmp_excess_quad <- tmp_excess


#### Simulate Excess Visits - Cubic -------------------------------------------

# pull out quadratic results
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
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_abx_patients)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>%
  mutate(excess_daysupp_pp = n_excess_abx_daysupp/n_excess_abx_patients,
         excess_daysupp_pp_all = n_excess_abx_daysupp/nrow(index_dates)) %>%  
  summarise(n_excess_patients_mean = mean(n_excess_abx_patients),
            n_excess_patients_median = median(n_excess_abx_patients),
            n_excess_patients_lo = quantile(n_excess_abx_patients, probs = 0.025),
            n_excess_patients_hi = quantile(n_excess_abx_patients, probs = 0.975),
            
            tot_excess_daysupp_mean = mean(n_excess_abx_daysupp),
            tot_excess_daysupp_median = median(n_excess_abx_daysupp),
            tot_excess_daysupp_lo = quantile(n_excess_abx_daysupp, probs = 0.025),
            tot_excess_daysupp_hi = quantile(n_excess_abx_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975))

# number of patients with excess abx
excess_patients <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure, 
         mean = n_excess_patients_mean,
         median = n_excess_patients_median,
         lo = n_excess_patients_lo,
         hi = n_excess_patients_hi) 

# fraction of patients with excess abx
frac_excess_patients <- excess_patients %>% 
  mutate(measure = "Percent of Patients with Excess ABX") %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates)))

# number of excess days supplied
excess_daysupp <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure, 
         mean = tot_excess_daysupp_mean,
         median = tot_excess_daysupp_median,
         lo = tot_excess_daysupp_lo,
         hi = tot_excess_daysupp_hi)

# number of excess days supplied per person
excess_daysupp_pp <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

# number of excess days supplied per person
excess_daysupp_pp_all <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure, 
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# compute excess numbers by antibiotic type
tmp <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$n_excess_patients_by_type)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  mutate(excess_daysupp_pp = n_excess_daysupp/n_excess_patients,
         excess_daysupp_pp_all = n_excess_daysupp/nrow(index_dates)) %>% 
  group_by(name) %>% 
  summarise(excess_patients_mean = mean(n_excess_patients),
            excess_patients_median = median(n_excess_patients),
            excess_patients_lo = quantile(n_excess_patients, probs = 0.025),
            excess_patients_hi = quantile(n_excess_patients, probs = 0.975),
            
            excess_daysupp_mean = mean(n_excess_daysupp),
            excess_daysupp_median = median(n_excess_daysupp),
            excess_daysupp_lo = quantile(n_excess_daysupp, probs = 0.025),
            excess_daysupp_hi = quantile(n_excess_daysupp, probs = 0.975),
            
            excess_daysupp_pp_mean = mean(excess_daysupp_pp),
            excess_daysupp_pp_median = median(excess_daysupp_pp),
            excess_daysupp_pp_lo = quantile(excess_daysupp_pp, probs = 0.025),
            excess_daysupp_pp_hi = quantile(excess_daysupp_pp, probs = 0.975),
            
            excess_daysupp_pp_all_mean = mean(excess_daysupp_pp_all),
            excess_daysupp_pp_all_median = median(excess_daysupp_pp_all),
            excess_daysupp_pp_all_lo = quantile(excess_daysupp_pp_all, probs = 0.025),
            excess_daysupp_pp_all_hi = quantile(excess_daysupp_pp_all, probs = 0.975)) 

# number of excess patients by type
excess_patients_by_type <- tmp %>% 
  mutate(measure = "Number of Patients with Excess ABX") %>% 
  select(measure,
         name,
         mean = excess_patients_mean,
         median = excess_patients_median,
         lo = excess_patients_lo,
         hi = excess_patients_hi)

# fraction of patients with excess abx by type
frac_excess_patients_by_type <- excess_patients_by_type %>% 
  mutate_at(vars(mean:hi),~100*(./nrow(index_dates))) %>% 
  mutate(measure = "Percent of Patients with Excess ABX")

# excess daysupplied by type
excess_daysupp_by_type <- tmp %>% 
  mutate(measure = "Total Excess Days Supplied") %>% 
  select(measure,
         name,
         mean = excess_daysupp_mean,
         median = excess_daysupp_median,
         lo = excess_daysupp_lo,
         hi = excess_daysupp_hi)

# excess days supplied per person by type
excess_daysupp_pp_by_type <- tmp %>% 
  mutate(measure = "Excess Days Supplied Per Patient (excess patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_mean,
         median = excess_daysupp_pp_median,
         lo = excess_daysupp_pp_lo,
         hi = excess_daysupp_pp_hi)

excess_daysupp_pp_all_by_type <- tmp %>%
  mutate(measure = "Excess Days Supplied Per Patient (all patients)") %>% 
  select(measure,
         name,
         mean = excess_daysupp_pp_all_mean,
         median = excess_daysupp_pp_all_median,
         lo = excess_daysupp_pp_all_lo,
         hi = excess_daysupp_pp_all_hi)

# bind_rows(excess_patients_by_type,
#           frac_excess_patients_by_type,
#           excess_daysupp_by_type,
#           excess_daysupp_pp_by_type,
#           excess_daysupp_pp_all_by_type) %>% 
#   mutate_at(vars(mean,lo,hi),~round(.,2)) %>% 
#   mutate(out = paste0(mean, " (",lo,", ",hi,")")) %>% 
#   select(measure,name,out) %>% 
#   spread(key = measure, value = out) %>% 
#   select(name,
#          `Number of Patients with Excess ABX`,
#          `Percent of Patients with Excess ABX`,
#          `Total Excess Days Supplied`,
#          `Excess Days Supplied Per Patient (excess patients)`,
#          `Excess Days Supplied Per Patient (all patients)`)

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
  summarise(mean = mean(n),
            median = median(n),
            lo = quantile(n,probs = 0.025),
            hi = quantile(n,probs = 0.975))


# mean number of excess antibiotics by patient across all patients
mean_excess_per_patient_all <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_patient_all)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (all patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

# mean number of excess antibiotics by patient across patients with excess
mean_excess_per_excess_patient <- tmp_excess %>% 
  mutate(res = map(excess_stats,~.$mean_excess_per_excess_patient)) %>% 
  select(trial,res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  summarise(measure = "Excess Number of ABX Prescriptions Per Patient (excess patients)",
            mean = mean(mean_excess),
            median = median(mean_excess),
            lo = quantile(mean_excess, probs = 0.025),
            hi = quantile(mean_excess, probs = 0.975))

## Aggregate stats
cubic_stats <- list(main_excess_stats = bind_rows(excess_patients,
                                                 frac_excess_patients,
                                                 excess_daysupp,
                                                 excess_daysupp_pp,
                                                 excess_daysupp_pp_all,
                                                 mean_excess_per_excess_patient,
                                                 mean_excess_per_patient_all),
                   
                   main_excess_stats_by_type = bind_rows(excess_patients_by_type,
                                                         frac_excess_patients_by_type,
                                                         excess_daysupp_by_type,
                                                         excess_daysupp_pp_by_type,
                                                         excess_daysupp_pp_all_by_type),
                   
                   excess_distribution = n_excess_abx_by_patient_count)

tmp_excess_cubic <- tmp_excess

#### Save Results ####

### Report info ----------------------------------------------------------------



save(index_dates,
     weekly_cp,
     final_abx_counts,
     abx_name_weekly_fits,
     p1,p2,p3,
     lm_res_plots,quad_res_plots,cubic_res_plots,
     p4,p5,
     excess_res1,excess_res2,
     lm_stats,quad_stats,cubic_stats,
     agg_boot_res,
     file = "~/Data/projects/excess_abx/histo/data/report_data.RData")

### Simulation results ---------------------------------------------------------

save(final_abx_counts,
     final_abx_visits,
     boot_res,
     tmp_excess_lm,
     tmp_excess_quad,
     tmp_excess_cubic,
     file = "~/Data/projects/excess_abx/histo/data/simulation_results.RData")

