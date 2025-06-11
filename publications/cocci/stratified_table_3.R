library(tidyverse)
library(lubridate)
library(smallDB)

# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "cocci"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 91 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

out_path <-paste0(delay_params$out_path,"risk_models/") 
sim_in_path <- paste0(delay_params$out_path,"sim_results/")
delay_base_path <- paste0(delay_params$base_path,"delay_results/")

# load AZ ind data
load(paste0(sim_in_path,"AZ_ind_data.RData"))

out_res <- tibble()

for (ind_x in c(3, 0, 1, NA)){
  
  if(!is.na(ind_x) && ind_x == 3){
    n_total <- nrow(AZ_ind_data)
    load(paste0(sim_in_path,"aggregated_sim_results.RData"))
  } else {
    load(paste0(sim_in_path,"aggregated_sim_results_AZ_", ind_x , ".RData"))
    if(is.na(ind_x)){
      n_total <- nrow(AZ_ind_data %>% filter(is.na(AZ)))
    } else {
      n_total <- nrow(AZ_ind_data %>% filter(AZ == ind_x))
    }
    
  }

  temp <- stringr::str_split(miss_bins_ssd$n, " ")
  n_new <- ceiling(as.numeric(sapply(temp, FUN = function(x) {x[1]})))
  n_range <- sapply(temp, FUN = function(x) {stringr::str_split(gsub("[()]", "", x[2]),  "-")})
  n_lower <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[1]})))
  n_upper <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[2]})))
  
  # No. of missed opportunities per patient
  
  miss_bins_ssd_new <- miss_bins_ssd %>% 
    mutate(miss_bin = paste0(">=",miss_bin)) %>% 
    select(miss_bin, n) %>% 
    mutate(N_pat =  n_new,
           N_pat_low = n_lower,
           N_pat_upper = n_upper) %>% 
    rowwise() %>% 
    mutate(estimate = paste0(format(N_pat, big.mark = ","), " (", format(round(N_pat/n_total *100, 1), nsmall = 1), "%)"),
           CI = paste0(format(N_pat_low, big.mark = ","), "-", format(N_pat_upper, big.mark = ","), " (",
                       format(round(N_pat_low/n_total *100, 1), nsmall = 1), "-",
                       format(round(N_pat_upper/n_total *100, 1), nsmall = 1), "%)")) %>% 
    select(metric = miss_bin, estimate, CI)
  
  miss_bins_ssd_new2 <- miss_bins_ssd %>% 
    mutate(miss_bin = paste0(">=",miss_bin)) %>% 
    select(miss_bin, n) %>% 
    mutate(N_pat =  n_new,
           N_pat_low = n_lower,
           N_pat_upper = n_upper) %>% 
    rowwise() %>% 
    mutate(estimate = paste0(format(round(N_pat/n_total *100, 1), nsmall = 1),  " (",
                             format(round(N_pat_low/n_total *100, 1), nsmall = 1), "-",
                             format(round(N_pat_upper/n_total *100, 1), nsmall = 1), ")")) %>% 
    select(metric = miss_bin, estimate)
  
  
  # mean/median No. of missed opportunities per patient
  miss_bins_ssds_mean_median <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_n_miss", "median_n_miss")) %>% 
    mutate(CI = paste0(low, "-", high),
           mean = as.character(mean)) %>% 
    select(metric = measure, estimate = mean, CI)
  
  miss_bins_ssds_mean_median2 <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_n_miss", "median_n_miss")) %>% 
    mutate(CI = paste0(low, "-", high),
           mean = as.character(mean)) %>% 
    mutate(estimate = paste0(mean, " (", CI,  ")") ) %>% 
    select(metric = measure, estimate)
  
  # Duration of Delayed Visits (days)
  
  temp <- stringr::str_split(dur_bins_ssd$n, " ")
  n_new <- ceiling(as.numeric(sapply(temp, FUN = function(x) {x[1]})))
  n_range <- sapply(temp, FUN = function(x) {stringr::str_split(gsub("[()]", "", x[2]),  "-")})
  n_lower <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[1]})))
  n_upper <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[2]})))
  
  dur_bins_ssd_new <- dur_bins_ssd %>%
    mutate(Duration = ifelse(duration_bin==1,
                             paste0(">= ",duration_bin, " Day"),
                             paste0(">= ",duration_bin, " Days"))) %>%
    select(Duration, n) %>%
    mutate(N_pat =  n_new,
           N_pat_low = n_lower,
           N_pat_upper = n_upper) %>%
    rowwise() %>%
    mutate(estimate = paste0(format(N_pat, big.mark = ","), " (", format(round(N_pat/n_total *100, 1), nsmall = 1), "%)"),
           CI = paste0(format(N_pat_low, big.mark = ","), "-", format(N_pat_upper, big.mark = ","), " (",
                       format(round(N_pat_low/n_total *100, 1), nsmall = 1), "-",
                       format(round(N_pat_upper/n_total *100, 1), nsmall = 1), "%)")) %>%
    select(metric = Duration, estimate, CI)
  
  dur_bins_ssd_new2 <- dur_bins_ssd %>% 
    mutate(Duration = ifelse(duration_bin==1,
                             paste0(">= ",duration_bin, " Day"),
                             paste0(">= ",duration_bin, " Days"))) %>% 
    select(Duration, n) %>% 
    mutate(N_pat =  n_new,
           N_pat_low = n_lower,
           N_pat_upper = n_upper) %>% 
    rowwise() %>% 
    mutate(estimate = paste0(format(round(N_pat/n_total *100, 1), nsmall = 1), " (",
                             format(round(N_pat_low/n_total *100, 1), nsmall = 1), "-",
                             format(round(N_pat_upper/n_total *100, 1), nsmall = 1), ")")) %>% 
    select(metric = Duration, estimate)
  
  
  # mean/median No. of missed opportunities per patient
  dur_bins_ssd_mean_median <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_dur", "median_dur")) %>% 
    mutate(CI = paste0(trimws(format(round(low, 2), nsmall = 2)), "-", trimws(format(round(high, 2), nsmall = 2))),
           mean = trimws(format(round(mean, 2), nsmall = 2))) %>% 
    select(metric = measure, estimate = mean, CI)
  
  dur_bins_ssd_mean_median2 <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_dur", "median_dur")) %>% 
    mutate(CI = paste0(trimws(format(round(low, 2), nsmall = 2)), "-", trimws(format(round(high, 2), nsmall = 2))),
           mean = trimws(format(round(mean, 2), nsmall = 2))) %>% 
    mutate(estimate = paste0(mean, " (", CI,  ")") ) %>% 
    select(metric = measure, estimate)
  
  
  delay_table <- bind_rows(tibble(metric = "No. of missed opportunities per patient",
                                  estimate = paste0("n total patients [AZ = ", ind_x, "] = ", n_total),
                                  CI = ""),
                           miss_bins_ssd_new %>% mutate(metric=paste0(" ", metric)),
                           miss_bins_ssds_mean_median %>% mutate(metric=paste0(" ", metric)),
                           tibble(metric = "Duration of Delayed Visits ",
                                  estimate = "",
                                  CI = ""),
                           dur_bins_ssd_new %>% mutate(metric=paste0(" ", metric)),
                           dur_bins_ssd_mean_median %>% mutate(metric=paste0(" ", metric)))
  
  delay_table2 <- bind_rows(tibble(metric = "No. of missed opportunities per patient",
                                   estimate = paste0("n total patients [AZ = ", ind_x, "] = ", n_total)),
                            miss_bins_ssd_new2 %>% mutate(metric=paste0(" ", metric)),
                            miss_bins_ssds_mean_median2 %>% mutate(metric=paste0(" ", metric)),
                            tibble(metric = "Duration of Delayed Visits ",
                                   estimate = ""),
                            dur_bins_ssd_new2 %>% mutate(metric=paste0(" ", metric)),
                            dur_bins_ssd_mean_median2 %>% mutate(metric=paste0(" ", metric)))
  
  
  if(!is.na(ind_x) && ind_x == 3) {
    names(delay_table2)[2] <- paste0("Estimate % of all patients [All patients]")
  }else {
    names(delay_table2)[2] <- paste0("Estimate % of all patients [AZ = ", ind_x, "]")
  }
  
      
  if(nrow(out_res) == 0){
    out_res <- delay_table2
  } else {
    out_res <- inner_join(out_res, delay_table2, by = "metric")
  }
  
  if(is.na(ind_x) || ind_x != 3) {
    write_csv(delay_table, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table3_AZ_", ind_x, ".csv"))
  }
}

write_csv(out_res, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table3_stratified.csv"))
