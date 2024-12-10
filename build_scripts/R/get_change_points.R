library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)
library(changepoint)
library(trend)
library(codeBuildr)



# load delay_parms
load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")

args = commandArgs(trailingOnly=TRUE)

# condition <- "dengue"
condition <- args[1]

delay_params <- delay_any_params[[condition]]


in_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",condition,"/delay_results")
# in_path <- paste0("~/Data/tmp/dengue/delay_results")

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",condition,"/change_point_results/")
# out_path <- "~/Data/Statepi_Diagnosis/prelim_results/dengue/change_point_results/"

if (!dir.exists(out_path)){
  dir.create(out_path)
}

# Load Data
load(paste0(in_path,"/all_dx_visits.RData"))

cond_name <- filter(codeBuildr::avail_disease_codes(),name==condition)$description
# cond_name <- filter(codeBuildr::avail_disease_codes(),name=="dengue")$description


##### Functions ####

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


#### Get Initial Counts ####


# count all visits
count_data_all <- visit_counts %>% 
  filter(is.na(dx_ver)) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7))

# load SSD codes
ssd_codes <- bind_rows(codeBuildr::load_ssd_codes(condition) %>%
                         filter(type == "icd9") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=9),
                       codeBuildr::load_ssd_codes(condition) %>%
                         filter(type == "icd10") %>%
                         select(dx=code) %>%
                         mutate(icd_ver=10))


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

#### Standard Approach - SSD Results -------------------------------------------

# piecewise model results
ssd_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_ssd,
                                cp_range = seq(from = 7, to = 300, by = 7),
                                model = .x,
                                periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))


# aggregate results output
ssd_res_out <- bind_rows(mutate(ssd_res1,
                                tmp = ifelse(periodicity," w/ periodicity",""),
                                label = paste0("Piecewise ", model, tmp))) %>% 
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

# split up results for plotting
tmp1 <- ssd_res1 %>%
  filter(periodicity)

tmp2 <- ssd_res1 %>%
  filter(!periodicity)

#### Optimal Change-point(s)

ssd_mod_rankings <- get_rankings(ssd_res_out) 

## Plots -----------------------------------------------------------------------

### Piecewise models with periodicity

plot_periodicity <- tmp1 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred2), color = "blue", size = .75) +
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp1, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp1, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()

### Piecewise models without periodicity

plot_no_periodicity <- tmp2 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred2), color = "blue", size = .75) +
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp2, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp2, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()


## Combine Results -------------------------------------------------------------

SSD_res_standard <- list(count_data_ssd = count_data_ssd,
                         ssd_res = ssd_res1,
                         ssd_res_out = ssd_res_out,
                         ssd_mod_rankings = ssd_mod_rankings,
                         plot_periodicity = plot_periodicity,
                         plot_no_periodicity = plot_no_periodicity)

rm(plot_periodicity, plot_no_periodicity,
   ssd_res1, ssd_res_out, ssd_mod_rankings,
   tmp1,tmp2, count_data_ssd)


#### Standard Approach - All Visits ####

# piecewise model results
all_res1 <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~find_cp(count_data_all,
                                cp_range = seq(from = 7, to = 300, by = 7),
                                model = .x,
                                periodicity = .y))) %>% 
  mutate(cp = map_dbl(cp_res,~.$cp)) %>% 
  mutate(pred = map(cp_res,~.$pred)) %>% 
  select(model,label,periodicity,cp,pred) %>% 
  mutate(pred_bound_cp = map_int(pred,find_pred_bound_cp)) %>% 
  mutate(label = paste0("Method: ",label,"\n Change point: ",cp," Pred Bound CP: ",pred_bound_cp))


# aggregate results output
all_res_out <- bind_rows(mutate(all_res1,
                                tmp = ifelse(periodicity," w/ periodicity",""),
                                label = paste0("Piecewise ", model, tmp))) %>% 
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

# split up results for plotting
tmp1 <- all_res1 %>%
  filter(periodicity)

tmp2 <- all_res1 %>%
  filter(!periodicity)

#### Optimal Change-point(s)

all_mod_rankings <- get_rankings(all_res_out) 

## Plots -----------------------------------------------------------------------

### Piecewise models with periodicity

plot_periodicity <- tmp1 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred2), color = "blue", size = .75) +
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp1, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp1, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()

### Piecewise models without periodicity

plot_no_periodicity <- tmp2 %>%
  select(label,cp,pred_bound_cp,pred) %>%
  unnest(pred) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse()+
  geom_line(aes(y = pred2), color = "blue", size = .75) +
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_vline(data = tmp2, aes(xintercept=cp), linetype=2) +
  geom_vline(data = tmp2, aes(xintercept=pred_bound_cp), linetype=2, color = "green") +
  facet_wrap(~label, ncol = 2) +
  theme_bw()


## Combine Results -------------------------------------------------------------

ALL_res_standard <- list(count_data_all = count_data_all,
                         all_res = all_res1,
                         all_res_out = all_res_out,
                         all_mod_rankings = all_mod_rankings,
                         plot_periodicity = plot_periodicity,
                         plot_no_periodicity = plot_no_periodicity)


rm(plot_periodicity, plot_no_periodicity,
   all_res1, all_res_out, all_mod_rankings,
   tmp1,tmp2, count_data_all)


### Save Results from Standard Approach ----------------------------------------


save(ALL_res_standard,SSD_res_standard,
     file = paste0(out_path,"standard_cp_res.RData"))




#### Bootstrapping Approach #####

# If change-point bounds are missing replace with default range from 10 to 100
if (is.na(delay_params$cp_lower) | is.na(delay_params$cp_upper)){
  cp_lower <- 10
  cp_upper <- 100
} else {
  cp_lower <- delay_params$cp_lower
  cp_upper <- delay_params$cp_upper
}

cp_range <- -seq(from = 7, to = cp_upper, by = 7)


ssd_set <- load_ssd_codes(condition) %>%
  mutate(dx_ver = ifelse(type == "icd9",9L,10L)) %>%
  select(dx = code, dx_ver)

upper_bound <- delay_params$upper_bound


###################
#### Prepars Data ####
###################

## Compute visit counts for entire population ----------------------------------

tmp1 <- all_dx_visits %>%
  inner_join(ssd_set,
             by = join_by(dx, dx_ver)) %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  mutate(code_set = "SSD")

tmp2 <- all_dx_visits %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  mutate(code_set = "All")

visit_counts <- bind_rows(tmp1,tmp2)
rm(tmp1,tmp2)

# add DOW factor
visit_counts <- visit_counts %>%
  filter(days_since_index<0) %>%
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow))


###################
#### Functions ####
###################

# This function computes bootstrapped visit counts, by generating a bootstrapped
# sample of patients and then computing their visit counts
compute_boot_count <- function(index_dates,dx_visit_data){
  sample_n(select(index_dates,patient_id,index_date), size = n(), replace = TRUE) %>%
    mutate(boot_id = row_number()) %>%
    inner_join(dx_visit_data,by = join_by(patient_id),relationship = "many-to-many") %>%
    filter(between(days_since_index,-upper_bound,1)) %>%
    distinct(boot_id,patient_id,days_since_index) %>%
    count(days_since_index) %>%
    filter(days_since_index<0) %>%
    mutate(dow = days_since_index %% 7) %>%
    mutate(dow = as.factor(dow))
}
# compute_boot_count(index_dx_dates,all_dx_visits)
# compute_boot_count(index_dx_dates,inner_join(all_dx_visits,ssd_set,by = join_by(dx, dx_ver)))


find_expected_fits <- function(data,cp) {
  
  if(cp>0){
    stop("CP must be <0")
  }
  
  fit_data <- dplyr::filter(data,days_since_index < cp)
  
  # fit trends for linear, quadratic, and cubic models -------------------------
  expected_fit_lm <- lm(n~days_since_index + dow, data = fit_data)
  expected_fit_quad <- lm(n~days_since_index+I(days_since_index^2) + dow, data = fit_data)
  expected_fit_cube <- lm(n~days_since_index+I(days_since_index^2)+I(days_since_index^3) + dow, data = fit_data)
  
  # Generate predictions for each model on whole data
  pred_fit_lm <- predict(expected_fit_lm, newdata = data, se = TRUE)
  pred_fit_quad <- predict(expected_fit_quad, newdata = data, se = TRUE)
  pred_fit_cube <- predict(expected_fit_cube, newdata = data, se = TRUE)
  
  # add fitted values to data along with SE for prediction bound ---------------
  expected_fit_data <- data %>%
    dplyr::mutate(linear = pred_fit_lm$fit) %>%
    dplyr::mutate(linear_se = sqrt(pred_fit_lm$se.fit^2 + pred_fit_lm$residual.scale^2)) %>%
    dplyr::mutate(quad = pred_fit_quad$fit) %>%
    dplyr::mutate(quad_se = sqrt(pred_fit_quad$se.fit^2 + pred_fit_quad$residual.scale^2)) %>%
    dplyr::mutate(cube = pred_fit_cube$fit) %>%
    dplyr::mutate(cube_se = sqrt(pred_fit_cube$se.fit^2 + pred_fit_cube$residual.scale^2)) %>%
    tidyr::gather(key = model, value = pred, -days_since_index, -n, -dow) %>%
    dplyr::mutate(pred_type = ifelse(str_detect(model,"se"),"pred_se","pred")) %>%
    dplyr::mutate(model = stringr::str_remove(model,"_se")) %>%
    tidyr::spread(key = pred_type, value = pred) %>%
    dplyr::arrange(model,days_since_index)
  
  out <- list(upper_bound = nrow(data),
              cp = cp,
              linear_model = coef(expected_fit_lm),
              quad_model = coef(expected_fit_quad),
              cube_model = coef(expected_fit_cube),
              expected_fits = expected_fit_data)
  
  return(out)
  
}

find_excess_fits <- function(expected_fit_data,cp){
  
  expected_fits <- expected_fit_data$expected_fits
  
  cp <- expected_fit_data$cp
  
  excess_fits <- expected_fits  %>%
    dplyr::mutate(excess = n-pred) %>%
    dplyr::filter(days_since_index >= cp) %>%
    dplyr::group_by(model) %>%
    tidyr::nest() %>%
    dplyr::mutate(fit = purrr::map(data,~loess(excess~days_since_index, data = .))) %>%
    dplyr::mutate(fit = purrr::map(fit,~predict(., se = TRUE))) %>%
    dplyr::mutate(data = purrr::map2(data,fit, ~dplyr::mutate(.x, excess_fit=.y$fit))) %>%
    dplyr::mutate(data = purrr::map2(data,fit, ~dplyr::mutate(.x, excess_se=sqrt(.y$se.fit^2 + .y$residual.scale^2)))) %>%
    dplyr::select(model,data) %>%
    tidyr::unnest(data) %>%
    dplyr::ungroup() %>%
    dplyr::select(model,days_since_index,excess:excess_se)
  
  expected_fit_data$excess_fits <- excess_fits
  
  return(expected_fit_data)
}

fit_cp <- function(data,cp){
  
  fit_data <- find_expected_fits(data,cp) %>%
    find_excess_fits()
  
  return(fit_data)
}

# This function fits the various trends for a given change-point to visit count data
assemble_cp_fit <- function(data,cp){
  
  fit_data <- fit_cp(data,cp)
  
  
  fit_data$expected_fits %>%
    select(model,days_since_index,dow,n,expected = pred) %>%
    left_join(select(fit_data$excess_fits,model,days_since_index,excess_fit),by = join_by(model, days_since_index)) %>%
    mutate(excess_fit = replace_na(excess_fit,0)) %>%
    rename(excess=excess_fit) %>%
    mutate(combined = expected+excess)
  
  
}

# This functions fits the various trends across a range of change-points
fit_trends <- function(count_data,cp_range = -7:-35){
  tibble(cp = cp_range) %>%
    mutate(fits = map(cp,~assemble_cp_fit(count_data, .))) %>%
    unnest(fits)
}


######################################
#### Fit Trends Across Bootstraps ####
######################################

## Generate Bootstrapped data --------------------------------------------------
boot_counts_all <- tibble(bootstrap = 1:100) %>%
  mutate(vis_counts = map(bootstrap,~compute_boot_count(index_dates = index_dx_dates,
                                                        dx_visit_data = all_dx_visits)))

boot_counts_ssd <- tibble(bootstrap = 1:100) %>%
  mutate(vis_counts = map(bootstrap,
                          ~compute_boot_count(index_dates = index_dx_dates,
                                              dx_visit_data = inner_join(all_dx_visits,
                                                                         ssd_set,
                                                                         by = join_by(dx, dx_ver)))))



## Fit Trends to Visit Count Data ----------------------------------------------

boot_fits_all <- boot_counts_all %>%
  mutate(fits = map(vis_counts,~fit_trends(.,cp_range = cp_range))) %>%
  select(bootstrap,fits)

boot_fits_ssd <- boot_counts_ssd %>%
  mutate(fits = map(vis_counts,~fit_trends(.,cp_range = cp_range))) %>%
  select(bootstrap,fits)

save(boot_fits_all,boot_fits_ssd, 
     file = paste0(out_path,"cp_boot_fits.RData"))


################################
#### Evaluate Change-Points ####
################################

## SSD Visits ------------------------------------------------------------------

# In-sample MSE across bootstraps
in_sample_mse_ssd <- boot_fits_ssd %>%
  unnest(fits) %>%
  group_by(bootstrap,cp,model) %>%
  summarise(mse = mean((combined-n)^2)) %>%
  group_by(cp,model) %>%
  summarise(in_mse = mean(mse))

# Out-of-sample MSE across bootstraps
out_of_sample_mse_ssd <- boot_fits_ssd %>%
  unnest(fits) %>%
  select(bootstrap,cp,model,days_since_index,combined) %>%
  inner_join(filter(visit_counts,code_set=="SSD")) %>%
  group_by(bootstrap,cp,model) %>%
  summarise(mse = mean((combined-n)^2)) %>%
  group_by(cp,model) %>%
  summarise(out_mse = mean(mse))

# Plot results
inner_join(in_sample_mse_ssd,out_of_sample_mse_ssd) %>%
  gather(key = key, value = value, -model, -cp) %>%
  ggplot(aes(cp,value,color = model)) +
  geom_line() +
  facet_wrap(~key,scale = "free_y")


cp_fits_ssd <- tibble(cp = cp_range) %>%
  mutate(fits = map(cp,~assemble_cp_fit(select(filter(visit_counts,code_set=="SSD"),-code_set),cp = .))) %>%
  unnest(fits)

## All Visits ------------------------------------------------------------------

# In-sample MSE across bootstraps
in_sample_mse_all <- boot_fits_all %>%
  unnest(fits) %>%
  group_by(bootstrap,cp,model) %>%
  summarise(mse = mean((combined-n)^2)) %>%
  group_by(cp,model) %>%
  summarise(in_mse = mean(mse))

# Out-of-sample MSE across bootstraps
out_of_sample_mse_all <- boot_fits_all %>%
  unnest(fits) %>%
  select(bootstrap,cp,model,days_since_index,combined) %>%
  inner_join(filter(visit_counts,code_set=="All")) %>%
  group_by(bootstrap,cp,model) %>%
  summarise(mse = mean((combined-n)^2)) %>%
  group_by(cp,model) %>%
  summarise(out_mse = mean(mse))

cp_fits_all <- tibble(cp = cp_range) %>%
  mutate(fits = map(cp,~assemble_cp_fit(select(filter(visit_counts,code_set=="All"),-code_set),cp = .))) %>%
  unnest(fits)

## Compute 100-fold boostrap CV of MSE across bootstraps -------------------------


## SSD visits
tmp <- boot_fits_ssd %>% 
  unnest(fits) %>% 
  select(bootstrap:days_since_index,n,combined)

boot_obs <- tmp %>% 
  distinct(bootstrap,days_since_index,n)


mse_res <- list()

for (i in 1:max(tmp$bootstrap)){
  
  tmp_mse <- tmp %>% 
    filter(bootstrap==i) %>% 
    select(days_since_index,cp,model,combined) %>% 
    inner_join(filter(boot_obs,bootstrap!=i),by = join_by(days_since_index),relationship = "many-to-many") %>% 
    group_by(bootstrap,model,cp) %>% 
    summarise(mse = mean((combined-n)^2))
  
  mse_res[[i]] <- tmp_mse
  
}

kfold_mse_ssd <- bind_rows(mse_res) %>% 
  group_by(model,cp) %>% 
  summarise(mean_mse = mean(mse),
            median_mse = median(mse)) %>% 
  ungroup()


## ALL visits
tmp <- boot_fits_all %>% 
  unnest(fits) %>% 
  select(bootstrap:days_since_index,n,combined)

boot_obs <- tmp %>% 
  distinct(bootstrap,days_since_index,n)


mse_res <- list()

for (i in 1:max(tmp$bootstrap)){
  
  tmp_mse <- tmp %>% 
    filter(bootstrap==i) %>% 
    select(days_since_index,cp,model,combined) %>% 
    inner_join(filter(boot_obs,bootstrap!=i),by = join_by(days_since_index),relationship = "many-to-many") %>% 
    group_by(bootstrap,model,cp) %>% 
    summarise(mse = mean((combined-n)^2))
  
  mse_res[[i]] <- tmp_mse
  
}

kfold_mse_all <- bind_rows(mse_res) %>% 
  group_by(model,cp) %>% 
  summarise(mean_mse = mean(mse),
            median_mse = median(mse)) %>% 
  ungroup()


### Save Output ----------------------------------------------------------------


save(out_of_sample_mse_ssd,out_of_sample_mse_all,
     in_sample_mse_ssd, in_sample_mse_all,
     kfold_mse_all,kfold_mse_ssd,
     cp_fits_all, cp_fits_ssd,
     visit_counts,
     file = paste0(out_path,"cp_boot_eval_res.RData"))


### Run Report -----------------------------------------------------------------



rmarkdown::render(input = "github/delay_diagnosis/build_scripts/R/report_scripts/opp_window_report.Rmd",
                  params = list(cond = condition),
                  output_dir = out_path)



