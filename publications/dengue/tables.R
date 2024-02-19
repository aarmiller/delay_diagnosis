
library(tidyverse)
library(lubridate)
library(smallDB)

# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "dengue"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
out_path <-paste0(delay_params$out_path,"risk_models/") 

# enrollment time prior to index
db <- src_sqlite(paste0(delay_params$small_db_path, cond_name, ".db"))

# Collect index dates
enroll_time_prior <- tbl(db,"index_dx_dates") %>% collect() %>% 
  distinct(patient_id, max_time_before_index) %>% 
  rename(enroll_time = max_time_before_index)

#load reg_demo data
load(paste0(out_path, "reg_data.RData"))

# Table 1 ----------------------------------------------------------------------
reg_demo <- reg_demo %>% inner_join(enroll_time_prior, by = "patient_id")


table1_fun <- function(demo_data){
  
  #age at dx
  
  age <- bind_rows(tibble("Characteristic" = "Age at diagnosis",
                      "Total Patients (% of patients)" = ""),
                      reg_demo %>% count(age_cat) %>% 
                     mutate(percent = n/(sum(n)) *100) %>% 
                     rowwise() %>% 
                     mutate(age_cat = paste0(" ", age_cat),
                            nice = paste0(format(n, big.mark = ","),
                                          " (", format(round(percent, 1), nsmall = 1),"%)")) %>% 
                     dplyr::select("Characteristic" = 1, "Total Patients (% of patients)" = nice),
                   tibble("Characteristic" = " Mean (SD)",
                          "Total Patients (% of patients)" = paste0(format(round(mean(reg_demo$age), 1), nsmall = 1), " (",
                                                                    format(round(sd(reg_demo$age), 1), nsmall = 1), ")")),
                   tibble("Characteristic" = " Median (IQR)",
                          "Total Patients (% of patients)" = paste0(format(round(median(reg_demo$age), 1), nsmall = 1), " (",
                                                                    format(round(IQR(reg_demo$age), 1), nsmall = 1), ")")),
  ) 
  
  #sex
  
  sex <- bind_rows(tibble("Characteristic" = "Sex",
                          "Total Patients (% of patients)" = ""),
                   reg_demo %>% count(female) %>%
                     mutate(female = factor(female, levels= c(F, T), label = c(" Male", " Female"))) %>% 
                     mutate(percent = n/(sum(n)) *100) %>% 
                     rowwise() %>% 
                     mutate(nice = paste0(format(n, big.mark = ","),
                                          " (", format(round(percent, 1), nsmall = 1),"%)")) %>% 
                     dplyr::select("Characteristic" = 1, "Total Patients (% of patients)" = nice)) 
  
  # data source
  
  data_source <- bind_rows(tibble("Characteristic" = "Database Source",
                                  "Total Patients (% of patients)" = ""),
                           reg_demo %>% count(source) %>%
                             mutate(source = factor(source, levels= c("ccae", "mdcr", "medicaid"),
                                                    label = c(" Commercial", " Supplemental Medicare", " Medicaid"))) %>% 
                             mutate(percent = n/(sum(n)) *100) %>% 
                             rowwise() %>% 
                             mutate(nice = paste0(format(n, big.mark = ","),
                                                  " (", format(round(percent, 1), nsmall = 1),"%)")) %>% 
                             dplyr::select("Characteristic" = 1, "Total Patients (% of patients)" = nice)) 
  
  # Enrollment time
  enroll_time <- bind_rows(tibble("Characteristic" = "Enrollment time prior to index (years)",
                                  "Total Patients (% of patients)" = ""),
                           tibble("Characteristic" = " Mean (SD)",
                                  "Total Patients (% of patients)" = paste0(format(round(mean(reg_demo$enroll_time/360), 1), nsmall = 1), " (",
                                                                            format(round(sd(reg_demo$enroll_time/360), 1), nsmall = 1), ")")),
                           tibble("Characteristic" = " Median (IQR)",
                                  "Total Patients (% of patients)" = paste0(format(round(median(reg_demo$enroll_time/360), 1), nsmall = 1), " (",
                                                                            format(round(IQR(reg_demo$enroll_time/360), 1), nsmall = 1), ")")),
                           tibble("Characteristic" = " Count <= 1 years",
                                  "Total Patients (% of patients)" = reg_demo %>% 
                                    filter(enroll_time <=360) %>% count() %>% 
                                    mutate(percent = n/(nrow(reg_demo)) *100) %>% 
                                    mutate(nice = paste0(format(n, big.mark = ","),
                                                         " (", format(round(percent, 1), nsmall = 1),"%)")) %>% .$nice),
                           # tibble("Characteristic" = " Count <= 1.5 years",
                           #        "Total Patients (% of patients)" = reg_demo %>% 
                           #          filter(enroll_time <=360+180) %>% count() %>% 
                           #          mutate(percent = n/(nrow(reg_demo)) *100) %>% 
                           #          mutate(nice = paste0(format(n, big.mark = ","),
                           #                               " (", format(round(percent, 1), nsmall = 1),"%)")) %>% .$nice),
                           tibble("Characteristic" = " Count <= 2 years",
                                  "Total Patients (% of patients)" = reg_demo %>% 
                                    filter(enroll_time <=360*2) %>% count() %>% 
                                    mutate(percent = n/(nrow(reg_demo)) *100) %>% 
                                    mutate(nice = paste0(format(n, big.mark = ","),
                                                         " (", format(round(percent, 1), nsmall = 1),"%)")) %>% .$nice),
                           tibble("Characteristic" = " Count <= 3 years",
                                  "Total Patients (% of patients)" = reg_demo %>% 
                                    filter(enroll_time <=360*3) %>% count() %>% 
                                    mutate(percent = n/(nrow(reg_demo)) *100) %>% 
                                    mutate(nice = paste0(format(n, big.mark = ","),
                                                         " (", format(round(percent, 1), nsmall = 1),"%)")) %>% .$nice),
                           tibble("Characteristic" = " Count > 3 years",
                                  "Total Patients (% of patients)" = reg_demo %>% 
                                    filter(enroll_time >360*3) %>% count() %>% 
                                    mutate(percent = n/(nrow(reg_demo)) *100) %>% 
                                    mutate(nice = paste0(format(n, big.mark = ","),
                                                         " (", format(round(percent, 1), nsmall = 1),"%)")) %>% .$nice)
  )
  
  bind_rows(age,sex,data_source, enroll_time)
}

table1 <- table1_fun(reg_demo)

write_csv(table1, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table1.csv"))


# Table 3 ----------------------------------------------------------------------
n_total <- nrow(reg_demo)
sim_in_path <- paste0(delay_params$out_path,"sim_results/")
load(paste0(sim_in_path,"aggregated_sim_results.RData"))

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

# mean/median No. of missed opportunities per patient
miss_bins_ssds_mean_median <- agg_stats_ssd$main_stats %>% 
  select(-measure_out) %>% 
  filter(measure %in% c("mean_n_miss", "median_n_miss")) %>% 
  mutate(CI = paste0(low, "-", high),
         mean = as.character(mean)) %>% 
  select(metric = measure, estimate = mean, CI)


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


# mean/median No. of missed opportunities per patient
dur_bins_ssd_mean_median <- agg_stats_ssd$main_stats %>% 
  select(-measure_out) %>% 
  filter(measure %in% c("mean_dur", "median_dur")) %>% 
  mutate(CI = paste0(low, "-", high),
         mean = as.character(mean)) %>% 
  select(metric = measure, estimate = mean, CI)

table3 <- bind_rows(tibble(metric = "No. of missed opportunities per patient",
                           estimate = "",
                           CI = ""),
                    miss_bins_ssd_new %>% mutate(metric=paste0(" ", metric)),
                    miss_bins_ssds_mean_median %>% mutate(metric=paste0(" ", metric)),
                    tibble(metric = "Duration of Delayed Visits ",
                           estimate = "",
                           CI = ""),
                    dur_bins_ssd_new %>% mutate(metric=paste0(" ", metric)),
                    dur_bins_ssd_mean_median %>% mutate(metric=paste0(" ", metric)))
                    
write_csv(table3, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table3.csv"))

                    
                    
