
library(tidyverse)
library(lubridate)
library(smallDB)

# devtools::install_github("aarmiller/smallDB")

# name of condition
proj_name <- "measles"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 14 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

out_path <-paste0(delay_params$out_path,"risk_models/") 
sim_in_path <- paste0(delay_params$out_path,"sim_results/")
delay_base_path <- paste0(delay_params$base_path,"delay_results/")

# enrollment time prior to index
db <- src_sqlite(paste0(delay_params$small_db_path, cond_name, ".db"))

# Collect index dates
enroll_time_prior <- tbl(db,"index_dx_dates") %>% collect() %>% 
  distinct(patient_id, max_time_before_index) %>% 
  rename(enroll_time = max_time_before_index)

#create reg_demo data
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),"index_cases.RData"))
load(paste0(delay_base_path,"demo_data.RData"))
demo1 <- demo1 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")
demo2 <- demo2 %>% select(-index_date) %>% 
  inner_join(index_cases %>% mutate(index_date = index_date + shift) %>% 
               select(patient_id, index_date), by = "patient_id")

reg_demo <- demo1 %>% 
  mutate(female=(sex==2),
         age = index_year-dobyr,
         stdrace = as.numeric(stdrace)) %>% 
  left_join(tibble(stdrace = c(0,1,2,4,9),
                   race = c("Missing/Unknown","White","Black","Hispanic","Other")),
            by = "stdrace") %>% 
  mutate(race = fct_relevel(race,"White"))

reg_demo <- reg_demo %>% 
  left_join(demo2 %>% 
              filter(index_date<=dtend & index_date>=dtstart) %>% 
              mutate(msa_new = msa %in% c("0","")) %>% 
              mutate(msa_new = ifelse(is.na(msa),NA,msa_new)) %>% 
              mutate(source = as.factor(source)) %>% 
              distinct(patient_id,source,msa=msa_new),
            by = "patient_id")

# age_cats <- c(-1,17,34,44,54,64,130)
age_cats <- c(-1,1,4,11,17,35,45,55,65,130) # pertussis age cats

reg_demo <- reg_demo %>% 
  mutate(age_cat = cut(age,age_cats)) %>% 
  mutate(month = as.character(month(as_date(index_date))),
         year = as.character(year(as_date(index_date)))) %>% 
  select(patient_id,female,age_cat,age,msa,source,year,month,msa,race) %>% 
  left_join(distinct(rural_visits,patient_id) %>% 
              mutate(rural = 1L), by = "patient_id") %>% 
  mutate(rural = replace_na(rural,0L))

reg_demo <- reg_demo %>% inner_join(enroll_time_prior, by = "patient_id")

# Table 1 ----------------------------------------------------------------------

table1_fun <- function(reg_demo){
  
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


# Table 2 ----------------------------------------------------------------------
# load(paste0(sim_in_path,"aggregated_sim_results.RData"))
# 
# index <- setting_counts_index_by_setting %>%
#   select(Setting,index_count_mean:index_pct2_upperCI) %>% 
#   mutate(across(index_count_mean:index_count_upperCI, ~formatC(ceiling(.), big.mark = ","))) %>% 
#   mutate(across(index_pct1_mean:index_pct2_upperCI, ~trimws(format(round(.,1), nsmall = 1)))) %>% 
#   mutate(index_count = paste0(index_count_mean," (",index_count_lowerCI,"-",index_count_upperCI,")"),
#          index_pct1 = paste0(index_pct1_mean," (",index_pct1_lowerCI,"-",index_pct1_upperCI,")")) %>% 
#   select(setting = Setting,`Index Visits`=index_count,
#          `% of all Index Locations`=index_pct1)
# 
# 
# misses <- setting_counts_index_by_setting %>%
#   select(Setting,n_mean:pct_opp_missed_high) %>% 
#   mutate(across(n_mean:n_high, ~formatC(ceiling(.), big.mark = ",")))%>% 
#   mutate(across(total_opps_mean:total_opps_high, ~formatC(ceiling(.), big.mark = ","))) %>% 
#   mutate(across(pct_opp_mean:pct_opp_high, ~trimws(format(round(.,1), nsmall = 1)))) %>% 
#   mutate(across(pct_opp_missed_mean:pct_opp_missed_high, ~trimws(format(round(.,1), nsmall = 1)))) %>% 
#   mutate(n = paste0(n_mean," (",n_low,"-",n_high,")"),
#          pct_opp = paste0(pct_opp_mean," (",pct_opp_low,"-",pct_opp_high,")"),
#          total_opp = paste0(total_opps_mean," (",total_opps_low,"-",total_opps_high,")"),
#          pct_opp_missed = paste0(pct_opp_missed_mean," (",pct_opp_missed_low,"-",pct_opp_missed_high,")"),) %>%
#   select(setting = Setting, 
#          `Missed Opportunities`=n,
#          `% Missed Opp. In Setting`=pct_opp,
#          'Total Diagnostic Opportunities' =total_opp,
#          `% of Opportunities Missed`=pct_opp_missed)
# 
# table2 <- inner_join(index, misses) %>% 
#   mutate(setting = factor(setting, levels = c("outpatient", "inpatient", "ed", "obs_stay"),
#                           labels = c("Outpatient", "Inpatient", "ED", "Observational Stay"))) %>% 
#   arrange(setting)

load(paste0(sim_in_path,"aggregated_sim_results.RData"))

table2 <- setting_table_new_raw_ssd %>% 
  mutate(`% of all Index Locations` = format(round(percent_index, 1), nsmall = 1)) %>% 
  mutate(`Number of Missed Opportunities` = paste0(miss_mean, "\n(", miss_low, "-", miss_high, ")")) %>% 
  mutate(`% of Missed opportunities` = paste0(format(round(frac_mean_all_miss_opp, 1), nsmall = 1), 
                                              "\n(", format(round(frac_low_all_miss_opp, 1), nsmall = 1), "-", format(round(frac_high_all_miss_opp, 1), nsmall = 1), ")")) %>% 
  mutate(`% of Missed Visit Days` = paste0(format(round(frac_mean, 1), nsmall = 1), "\n(", format(round(frac_low, 1), nsmall = 1), "-", 
                                           format(round(frac_high, 1), nsmall = 1), ")")) %>% 
  select(Setting,
         `Index Count`=index_n, 
         `% of all Index Locations`,
         `Potential Missed Opportunity Visit Days`=potential_opps,
         `Number of Missed Opportunities`, 
         `% of Missed opportunities`,
         `% of Missed Visit Days`) %>% 
  mutate(`% of Missed Visit Days` = ifelse(Setting == "inpatient visit", "", `% of Missed Visit Days`),
         `% of Missed opportunities` = ifelse(Setting == "inpatient visit", "", `% of Missed opportunities`)) %>% 
  select(1,4,2,3,5,6) %>% 
  mutate(across(where(is.integer), ~ formatC(.x,  big.mark = ",")))
write_csv(table2, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/table2.csv"))



# Table 3 ----------------------------------------------------------------------
n_total <- nrow(reg_demo)

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
  mutate(CI = paste0(trimws(format(round(low, 2), nsmall = 2)), "-", trimws(format(round(high, 2), nsmall = 2))),
         mean = trimws(format(round(mean, 2), nsmall = 2))) %>% 
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

# Appendix Table 1 --------------------------------------------------------------
# ICD code SSDs
load("/Volumes/AML/params/final_delay_params.RData")
proj_name <- "measles"
delay_params <- final_delay_params[[proj_name]]

codeBuildr:::load_ssd_codes(delay_params$ssd_name) %>% 
  mutate(icd_version = ifelse(type == "icd9", 9L, 10L)) %>% 
  select(icd_codes = code, icd_version) %>% 
  inner_join(codeBuildr::labels) %>% 
  mutate(icd_codes_nice = map_chr(icd_codes, ~icd::short_to_decimal(.))) %>% 
  mutate(icd_codes_nice = paste0("'",icd_codes_nice ))%>% 
  select(icd_codes, icd_codes_nice, icd_version, desc) %>% 
write_csv(paste0("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/appendix_table1.csv"))

# Appendix Table 3 --------------------------------------------------------------
# sensitivity analysis all visits

n_total <- nrow(reg_demo)

temp <- stringr::str_split(miss_bins_all$n, " ")
n_new <- ceiling(as.numeric(sapply(temp, FUN = function(x) {x[1]})))
n_range <- sapply(temp, FUN = function(x) {stringr::str_split(gsub("[()]", "", x[2]),  "-")})
n_lower <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[1]})))
n_upper <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[2]})))

# No. of missed opportunities per patient

miss_bins_all_new <- miss_bins_all %>% 
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
miss_bins_alls_mean_median <- agg_stats_all$main_stats %>% 
  select(-measure_out) %>% 
  filter(measure %in% c("mean_n_miss", "median_n_miss")) %>% 
  mutate(CI = paste0(low, "-", high),
         mean = as.character(mean)) %>% 
  select(metric = measure, estimate = mean, CI)


# Duration of Delayed Visits (days)

temp <- stringr::str_split(dur_bins_all$n, " ")
n_new <- ceiling(as.numeric(sapply(temp, FUN = function(x) {x[1]})))
n_range <- sapply(temp, FUN = function(x) {stringr::str_split(gsub("[()]", "", x[2]),  "-")})
n_lower <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[1]})))
n_upper <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[2]})))


dur_bins_all_new <- dur_bins_all %>% 
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
dur_bins_all_mean_median <- agg_stats_all$main_stats %>% 
  select(-measure_out) %>% 
  filter(measure %in% c("mean_dur", "median_dur")) %>% 
  mutate(CI = paste0(trimws(format(round(low, 2), nsmall = 2)), "-", trimws(format(round(high, 2), nsmall = 2))),
         mean = trimws(format(round(mean, 2), nsmall = 2))) %>% 
  select(metric = measure, estimate = mean, CI)

table3 <- bind_rows(tibble(metric = "No. of missed opportunities per patient",
                           estimate = "",
                           CI = ""),
                    miss_bins_all_new %>% mutate(metric=paste0(" ", metric)),
                    miss_bins_alls_mean_median %>% mutate(metric=paste0(" ", metric)),
                    tibble(metric = "Duration of Delayed Visits ",
                           estimate = "",
                           CI = ""),
                    dur_bins_all_new %>% mutate(metric=paste0(" ", metric)),
                    dur_bins_all_mean_median %>% mutate(metric=paste0(" ", metric)))

write_csv(table3, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/appendix_table3.csv"))


# Appendix Table 2 --------------------------------------------------------------
# sensitivity analysis varying delay window

appendix_table2 <- tibble()

for(i in final_delay_params[[proj_name]]$cp){
  delay_params$cp <- i
  delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")
  
  sim_in_path <- paste0(delay_params$out_path,"sim_results/")
  load(paste0(sim_in_path,"aggregated_sim_results.RData"))
  
  temp <- stringr::str_split(miss_bins_ssd$n, " ")
  n_new <- ceiling(as.numeric(sapply(temp, FUN = function(x) {x[1]})))
  n_range <- sapply(temp, FUN = function(x) {stringr::str_split(gsub("[()]", "", x[2]),  "-")})
  n_lower <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[1]})))
  n_upper <-  ceiling(as.numeric(sapply(n_range, FUN = function(x) {x[2]})))
  
  percent_pateint_delayed <- miss_bins_ssd %>% 
    mutate(miss_bin = paste0(">=",miss_bin)) %>% 
    select(miss_bin, n) %>% 
    mutate(N_pat =  n_new,
           N_pat_low = n_lower,
           N_pat_upper = n_upper) %>% 
    rowwise() %>% 
    mutate(nice = paste0(format(round(N_pat/n_total *100, 1), nsmall = 1),
                         "\n(95%CI: ",  format(round(N_pat_low/n_total *100, 1), nsmall = 1),
                         "-",
                         format(round(N_pat_upper/n_total *100, 1), nsmall = 1),")")) %>% 
    filter(miss_bin == ">=1") %>% 
    mutate(miss_bin = "percent_pateint_delayed") %>% 
    select(metric = miss_bin, nice)
  
  names(percent_pateint_delayed)[2] <- paste0("delay_window_1_", i-1)
  
  # mean/median No. of missed opportunities per patient
  miss_bins_ssds_mean_median <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_n_miss")) %>% 
    mutate(CI = paste0(low, "-", high),
           mean = as.character(mean)) %>% 
    mutate(nice = paste0(mean, "\n(95%CI: ", CI, ")")) %>% 
    select(metric = measure, nice) 
  
  names(miss_bins_ssds_mean_median)[2] <- paste0("delay_window_1_", i-1)
  
  # mean/median No. of missed opportunities per patient
  dur_bins_ssd_mean_median <- agg_stats_ssd$main_stats %>% 
    select(-measure_out) %>% 
    filter(measure %in% c("mean_dur")) %>% 
    mutate(CI = paste0(trimws(format(round(low, 2), nsmall = 2)), "-", trimws(format(round(high, 2), nsmall = 2))),
           mean = trimws(format(round(mean, 2), nsmall = 2))) %>% 
    mutate(nice = paste0(mean, "\n(95%CI: ", CI, ")")) %>% 
    select(metric = measure, nice) 
  
  names(dur_bins_ssd_mean_median)[2] <- paste0("delay_window_1_", i-1)
  
  temp_table <- bind_rows(percent_pateint_delayed,
                          miss_bins_ssds_mean_median,
                          dur_bins_ssd_mean_median)
  if(nrow(appendix_table2) == 0){
    appendix_table2 <- temp_table
  } else {
    appendix_table2 <- inner_join(appendix_table2, temp_table)
  }
  
}

write_csv(appendix_table2, paste0("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", proj_name, "/tables/appendix_table2.csv"))
