

rm(list = ls())
library(tidyverse)
library(smallDB)
library(codeBuildr)

########################
#### Initial Params ####
########################

args = commandArgs(trailingOnly=TRUE)

# condition <- "cocci"
condition <- args[1]

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",condition,"/change_point_results/")

load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")
delay_params <- delay_any_params[[condition]]

# If change-point bounds are missing replace with default range from 10 to 100
if (is.na(delay_params$cp_lower) | is.na(delay_params$cp_upper)){
  cp_lower <- 10
  cp_upper <- 100
} else {
  cp_lower <- delay_params$cp_lower
  cp_upper <- delay_params$cp_upper
}

cp_range <- (-cp_lower):(-cp_upper)

ssd_set <- load_ssd_codes(condition) %>%
  mutate(dx_ver = ifelse(type == "icd9",9L,10L)) %>%
  select(dx = code, dx_ver)

upper_bound <- delay_params$upper_bound


###################
#### Load Data ####
###################

### Collect data to generate initial trends ------------------------------------

# Here we will start by collecting the basic counts of all visits before dengue
# diagnosis

load(paste0("/Shared/Statepi_Diagnosis/prelim_results/",condition,"/delay_results/all_dx_visits.RData"))
# load(paste0("/Volumes/Statepi_Diagnosis/prelim_results/",condition,"/delay_results/all_dx_visits.RData"))

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

# Compare Boostrapped data to real data
# boot_counts_all %>%
#   unnest(vis_counts) %>%
#   left_join(filter(visit_counts,code_set=="All") %>%
#               mutate(bootstrap=1) %>%
#               select(bootstrap,days_since_index,n_true = n)) %>%
#   ggplot(aes(days_since_index,n,group = bootstrap)) +
#   geom_line() +
#   geom_line(aes(y = n_true),color = "red", size = 0.7) +
#   theme_bw()
#
# boot_counts_ssd %>%
#   unnest(vis_counts) %>%
#   left_join(filter(visit_counts,code_set=="SSD") %>%
#               mutate(bootstrap=1) %>%
#               select(bootstrap,days_since_index,n_true = n)) %>%
#   ggplot(aes(days_since_index,n,group = bootstrap)) +
#   geom_line() +
#   geom_line(aes(y = n_true),color = "red", size = 0.7) +
#   theme_bw()

## Fit Trends to Visit Count Data ----------------------------------------------

boot_fits_all <- boot_counts_all %>%
  mutate(fits = map(vis_counts,~fit_trends(.,cp_range = cp_range))) %>%
  select(bootstrap,fits)

# boot_fits_all$fits[[1]] %>%
#   ggplot(aes(days_since_index,n)) +
#   geom_line(aes(y = combined, group = cp)) +
#   geom_point()

boot_fits_ssd <- boot_counts_ssd %>%
  mutate(fits = map(vis_counts,~fit_trends(.,cp_range = cp_range))) %>%
  select(bootstrap,fits)

# boot_fits_ssd$fits[[1]] %>%
#   ggplot(aes(days_since_index,n)) +
#   geom_line(aes(y = combined, group = cp)) +
#   geom_point()

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

# Find Best Model
# out_of_sample_mse_ssd %>%
#   ungroup() %>%
#   arrange(out_mse)

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

# Plot results
# inner_join(in_sample_mse_all,out_of_sample_mse_all) %>%
#   gather(key = key, value = value, -model, -cp) %>%
#   ggplot(aes(cp,value,color = model)) +
#   geom_line() +
#   facet_wrap(~key,scale = "free_y")

# Find Best Model
# out_of_sample_mse_all %>%
#   ungroup() %>%
#   arrange(out_mse)

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
  
  print(i)
  
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
  
  print(i)
  
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



