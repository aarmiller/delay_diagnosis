


### Extract All Visits ###

## Run on Argon ----------------------------------------------------------------

library(tidyverse)
library(icd)

pd_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/pd/pd.db")
als_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/als/als.db")
mg_db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/mg/mg.db")

pd_dx_vis <- pd_db %>% 
  tbl("all_dx_visits") %>% 
  filter(days_since_index<=0) %>% 
  collect() %>% 
  distinct(patient_id,dx,dx_ver,days_since_index)

als_dx_vis <- als_db %>% 
  tbl("all_dx_visits") %>% 
  filter(days_since_index<=0) %>% 
  collect() %>% 
  distinct(patient_id,dx,dx_ver,days_since_index)

mg_dx_vis <- mg_db %>% 
  tbl("all_dx_visits") %>% 
  filter(days_since_index<=0) %>% 
  collect() %>% 
  distinct(patient_id,dx,dx_ver,days_since_index)


pd_index <- pd_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

als_index <- als_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

mg_index <- mg_db %>% 
  tbl("index_dx_dates") %>% 
  collect()

all_dx_vis <- bind_rows(mutate(pd_dx_vis,group = "PD"),
                         mutate(als_dx_vis,group = "ALS"),
                         mutate(mg_dx_vis,group = "MG"))

all_index <- bind_rows(mutate(pd_index,group = "PD"),
                       mutate(als_index,group = "ALS"),
                       mutate(mg_index,group = "MG"))

save(all_dx_vis,all_index,
     file = "/Shared/AML/tmp_transfer/jacob_pd_grant_all_dx_index.RData")

#####################
#### Prepare Data ###
#####################

# Run Locally

rm(list = ls())

## Compute SSD Counts ----------------------------------------------------------

load("/Volumes/AML/tmp_transfer/jacob_pd_grant_all_dx_index.RData")

# Load SSD Codes
ssd_codes <- read_csv("~/Downloads/top_500.csv")

## Clean Codes -----------------------------------------------------------------
ssd_codes <- ssd_codes %>% 
  select(icd_9, coding) %>% 
  mutate(icd_9 = stringr::str_replace(icd_9, "\\.", "")) %>% 
  mutate(icd_9 = case_when(icd::is_billable(icd_9) ~ icd_9,
                           !icd::is_billable(icd_9) ~ glue::glue("0{icd_9}")
                           )
         ) 

ssd_codes <- ssd_codes %>% 
  rename(dx = icd_9) %>% 
  mutate(dx_ver = 9L)

# SSD visits
all_ssd_vis <- all_dx_vis %>% 
  inner_join(select(ssd_codes,dx,dx_ver)) %>% 
  distinct(group,patient_id,dx,days_since_index)

# Individual Counts
ssd_counts <- all_ssd_vis %>% 
  group_by(group,dx,days_since_index) %>% 
  count(dx)

# Any Count
any_ssd_counts <- all_ssd_vis %>% 
  distinct(group,patient_id,days_since_index) %>% 
  count(group,days_since_index)


############################
#### Analysis Functions ####
############################

find_expected_fits <- function(data,cp) {
  
  if(cp>0){
    stop("CP must be <0")
  }
  
  fit_data <- dplyr::filter(data,days_since_index < cp)
  
  # fit trends for linear, quadratic, and cubic models -------------------------
  expected_fit_lm <- lm(n~days_since_index, data = fit_data)
  expected_fit_quad <- lm(n~days_since_index+I(days_since_index^2), data = fit_data)
  expected_fit_cube <- lm(n~days_since_index+I(days_since_index^2)+I(days_since_index^3), data = fit_data)
  
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
    tidyr::gather(key = model, value = pred, -days_since_index, -n) %>%
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
    select(model,days_since_index,n,expected = pred) %>%
    left_join(select(fit_data$excess_fits,model,days_since_index,excess_fit),by = join_by(model, days_since_index)) %>%
    mutate(excess_fit = replace_na(excess_fit,0)) %>%
    rename(excess=excess_fit) %>%
    mutate(combined = expected+excess)
  
  
}


######################
#### Analyze Data ####
######################


cp <- -160
obs_window <- 4*365

all_index %>% 
  filter(group == "PD") %>% 
  filter(time_before_index>=obs_window) %>% 
  distinct(patient_id)



count_data <- all_ssd_vis %>% 
  inner_join(all_index %>% 
               filter(group == "PD") %>% 
               filter(time_before_index>=obs_window) %>% 
               distinct(patient_id)) %>% 
  distinct(group,patient_id,days_since_index) %>% 
  count(group,days_since_index) %>% 
  filter(between(days_since_index,-obs_window,-1)) %>% 
  filter(group == "PD") %>% 
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(days_since_index,n,dow)

count_data %>% 
  mutate(period = days_since_index %/% 7) %>% 
  group_by(period) %>% 
  summarise(n = sum(n)) %>% 
  ggplot(aes(period,n)) +
  geom_point()


tmp <- count_data %>% 
  arrange(days_since_index) %>% 
  mutate(ma_n = (n + lag(n) + lag(n,2) +
              lag(n,3) + lag(n,4) +
              lag(n,5) + lag(n,6))/7)

tmp %>% 
  ggplot(aes(days_since_index,ma_n)) +
  geom_point()

count_data %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_point() +
  theme_bw() +
  scale_x_reverse() +
  ggtitle("All SSD visits prior to PD") +
  xlab("Days Before Diagnosis") +
  ylab("Visit Count")

count_data %>% 
  mutate(period = days_since_index %/% 7) %>% 
  group_by(period) %>% 
  summarise(n = sum(n)) %>% 
  # filter(period>-157) %>% 
  mutate(week = -period) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  theme_bw() +
  scale_x_reverse() +
  ggtitle("Weekly Total of All SSD visits prior to PD") +
  xlab("Weeks Before Diagnosis") +
  ylab("Visit Count")

tmp %>% 
  ggplot(aes(-days_since_index,ma_n)) +
  geom_point() + 
  scale_x_reverse() +
  theme_bw() +
  ggtitle("Moving Average of All SSD visits prior to PD") +
  xlab("Days Before Diagnosis") +
  ylab("Visit Count")


any_ssd_counts %>% 
  filter(obs_window == 365*3) %>% 
  filter(group == "PD") %>% 
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(days_since_index,n,dow) %>% 
  assemble_cp_fit(cp = cp) %>% 
  filter(model == "quad") %>% 
  mutate(combined = ifelse(days_since_index > cp,combined,NA)) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_point() +
  geom_line(aes(y=expected), color = "blue", size = 1.5) +
  theme_bw() +
  geom_line(aes(y = combined), color = "red", size = 1.5) +
  geom_ribbon(aes(ymin = expected, ymax = combined), fill = "red", alpha = 0.5)


### Final Count Data -------------------

obs_window <- 4*365

count_data <- all_ssd_vis %>% 
  inner_join(all_index %>% 
               filter(group == "PD") %>% 
               filter(time_before_index>=obs_window) %>% 
               distinct(patient_id)) %>% 
  distinct(group,patient_id,days_since_index) %>% 
  count(group,days_since_index) %>% 
  filter(between(days_since_index,-obs_window,-1)) %>% 
  filter(group == "PD") %>% 
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(days_since_index,n,dow)

final_count_data <- count_data %>% 
  mutate(period = days_since_index %/% 7) %>% 
  group_by(period) %>% 
  summarise(n = sum(n))

final_count_data <- final_count_data %>% 
  filter(period > -209)

tmp <- final_count_data %>% 
  rename(days_since_index = period)

cp <- -43

pop_size <- all_index %>% 
  filter(group == "PD") %>% 
  filter(time_before_index>=obs_window) %>% 
  distinct(patient_id) %>% 
  nrow()

plot_data <- assemble_cp_fit(tmp,cp = cp) %>% 
  mutate(combined = ifelse(days_since_index > cp,combined,NA)) %>% 
  mutate_at(vars(n:combined), ~1000*(./pop_size))

save(plot_data, file = "~/OneDrive - University of Iowa/grant_proposals/jacob_pd/data/main_plot_data.RData")


plot_data %>% 
  filter(model=="quad") %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_point(size = 2.5) +
  geom_line(aes(y=expected), color = "blue", size = 1.5, linetype = 1) +
  theme_bw() +
  geom_line(aes(y = combined), color = "red", size = 1.5) +
  geom_ribbon(aes(ymin = expected, ymax = combined), fill = "red", alpha = 0.5) +
  scale_x_reverse(breaks=c(52,104,156)) +
  geom_vline(aes(xintercept = -cp), linetype = 2) +
  annotate("text", x = 20, y = 580, label = "Diagnostic \n Opportunity Window",
           fontface = 2, size = 6) +
  coord_cartesian(ylim = c(240,540), clip = "off") +
  theme(plot.margin = margin(55, 40, 10, 10)) +
  xlab("Weeks Before Diagnosis") +
  ylab("Number of SSD Visits (per 1K patients)") +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))
ggsave("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/figures/main_delay_plot.svg", 
       height = 6,width = 8)
ggsave("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/figures/main_delay_plot.pdf", 
       height = 6,width = 8)



?annotate





# This function fits the various trends for a given change-point to visit count data
assemble_cp_fit <- function(data,cp){
  
  fit_data <- fit_cp(data,cp)
  
  
  fit_data$expected_fits %>%
    select(model,period,dow,n,expected = pred) %>%
    left_join(select(fit_data$excess_fits,model,period,excess_fit),by = join_by(model, period)) %>%
    mutate(excess_fit = replace_na(excess_fit,0)) %>%
    rename(excess=excess_fit) %>%
    mutate(combined = expected+excess)
  
  
}


cp <- -160

any_ssd_counts %>% 
  filter(obs_window == 365*3) %>% 
  filter(group == "PD") %>% 
  mutate(dow = period %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(period,n,dow) %>% 
  assemble_cp_fit(cp = cp) %>% 
  filter(model == "quad") %>% 
  mutate(combined = ifelse(period > cp,combined,NA)) %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y=expected), color = "blue", size = 1.5) +
  theme_bw() +
  geom_line(aes(y = combined), color = "red", size = 1.5) +
  geom_ribbon(aes(ymin = expected, ymax = combined), fill = "red", alpha = 0.5)
