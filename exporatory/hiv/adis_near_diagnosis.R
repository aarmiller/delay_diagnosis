

rm(list = ls())

library(tidyverse)
library(codeBuildr)



###################
#### Functions ####
###################

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


###################
#### Code Sets ####
###################


#### ADI - Aids-defining illnesses ---------------------------------------------

adi_codes <- list(
                  # oral_hair_leycoplakia = list(icd9 = c("5286"),
                  #                              icd10 = c("K133")),
                  
                  candidiasis_of_esophagus = list(icd9 = c("1128"),
                                                  icd10 = c("C53")),
                  
                  cryptococcosis = list(icd9 = c("321"),
                                        icd10 = c("A072")),
                  
                  cryptosporidiosis = list(icd9 = c("0074"),
                                           icd10 = c("B25")),
                  
                  # cytomegalovirus = list(icd9 = c("4841"),
                  #                        icd10 = c("G934")),
                  
                  encephalopathy = list(icd9 = c("3483"),
                                        icd10 = c("B002")),
                  
                  # histoplasmosis = list(icd9 = c("115"),
                  #                       icd10 = c("A073")),
                  
                  # isosporiasis = list(icd9 = c("0072"),
                  #                     icd10 = c("A073")),
                  
                  kaposi_sarcoma = list(icd9 = c("176"),
                                        icd10 = c("C46")),
                  
                  pneumocystosis = list(icd9 = c("1363"),
                                        icd10 = c("B59")),
                  
                  leukoencephalopathy = list(icd9 = c("463"),
                                             icd10 = c("A812")),
                  
                  meningoencephalitis_toxo = list(icd9 = c("130"),
                                                  icd10 = c("B582")))

tmp1 <- enframe(map(adi_codes,~.$icd9),name = "label",value = "dx") %>% 
  unnest(dx)

tmp1 <- tmp1 %>% 
  mutate(dx_expand = map(dx,children_safe)) %>% 
  unnest(dx_expand) %>% 
  distinct(label,dx = dx_expand) %>% 
  mutate(dx_ver = 9L)


tmp2 <- enframe(map(adi_codes,~.$icd10),name = "label",value = "dx") %>% 
  unnest(dx)

tmp2 <- tmp2 %>% 
  mutate(dx_expand = map(dx,children_safe)) %>% 
  unnest(dx_expand) %>% 
  distinct(label,dx = dx_expand) %>% 
  mutate(dx_ver = 10L)

final_adi_codes <- bind_rows(tmp1,tmp2) %>% 
  rename(phil_label = label) %>% 
  inner_join(codeBuildr::all_icd_labels)

# write_csv(final_adi_codes,"~/Downloads/ADI_codes.csv")

#### SSD codes -----------------------------------------------------------------

ssd_codes <- list(Diarrhea = list(icd9 = c("78791"),
                                  icd10 = c("R197")),
                  
                  `Thrombocytopenia / Anemia` = list(icd9 = c("2875","2859"),
                                                      icd10 = c("D696","D649")),
                  
                  Pneumonia = list(icd9 = c("481","486"),
                                                icd10 = c("J189","J13")),
                  
                  `Enlarged Lymph Nodes` = list(icd9 = c("7856"),
                                     icd10 = c("R591")),
                  
                  Weightloss = list(icd9 = c("78321"),
                                    icd10 = c("R634")),
                  
                  Rash = list(icd9 = c("7821"),
                              icd10 = c("R21"))
                  )

tmp1 <- enframe(map(ssd_codes,~.$icd9),name = "label",value = "dx") %>% 
  unnest(dx) 

tmp1 <- tmp1 %>% 
  mutate(dx_expand = map(dx,children_safe)) %>% 
  unnest(dx_expand) %>% 
  distinct(label,dx = dx_expand) %>% 
  mutate(dx_ver = 9L)


tmp2 <- enframe(map(ssd_codes,~.$icd10),name = "label",value = "dx") %>% 
  unnest(dx)

tmp2 <- tmp2 %>% 
  mutate(dx_expand = map(dx,children_safe)) %>% 
  unnest(dx_expand) %>% 
  distinct(label,dx = dx_expand) %>% 
  mutate(dx_ver = 10L)

final_ssd_codes <- bind_rows(tmp1,tmp2) 


##################
#### Analysis ####
##################

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/hiv/hiv.db")

index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()

all_dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  collect()

adi_visits <- all_dx_visits %>% 
  inner_join(distinct(final_adi_codes,dx,dx_ver)) 

ssd_visits <- all_dx_visits %>% 
  inner_join(distinct(final_ssd_codes,dx,dx_ver,label)) %>% 
  distinct(patient_id,days_since_index,label)

# adi_window <- 30
# 
# adi_ids <- adi_visits %>% 
#   filter(between(days_since_index,-adi_window,adi_window)) %>% 
#   distinct(patient_id)
# 
# ssd_counts <- ssd_visits %>% 
#   filter(between(days_since_index,-365*2,-1)) %>% 
#   left_join(mutate(adi_ids,adi_patient=TRUE)) %>% 
#   mutate(adi_patient=replace_na(adi_patient,FALSE)) %>% 
#   group_by(adi_patient) %>% 
#   count(days_since_index) %>% 
#   ungroup() 
# 
# ssd_counts %>% 
#   ggplot(aes(days_since_index,n)) +
#   geom_line() +
#   facet_wrap(~adi_patient,scales = "free_y")


### Fit SSD trends -------------------------------------------------------------

tmp_adi <- filter(ssd_counts,adi_patient) %>%
  filter(days_since_index<0) %>%
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(-adi_patient)

tmp_non_adi <- filter(ssd_counts,!adi_patient) %>%
  filter(days_since_index<0) %>%
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow)) %>% 
  select(-adi_patient)

adi_fits <- fit_trends(tmp_adi,cp_range = -250:(-100))
non_adi_fits <- fit_trends(tmp_non_adi,cp_range = -250:(-100))

adi_fits %>% 
  # filter(days_since_index<cp) %>% 
  group_by(cp,model) %>% 
  summarise(mse = mean((n-combined)^2)) %>% 
  ungroup() %>% 
  filter(model == "linear") %>%
  arrange(mse)

non_adi_fits %>% 
  # filter(days_since_index<cp) %>% 
  group_by(cp,model) %>% 
  summarise(mse = mean((n-combined)^2)) %>% 
  ungroup() %>% 
  filter(model == "linear") %>%
  arrange(mse)


### Prelim Counts --------------------------------------------------------------

index_dx_dates %>% 
  summarise(n_1yr = sum(time_before_index>=365),
            n_2yr = sum(time_before_index>=365*2),
            n_3yr = sum(time_before_index>=365*3),
            n_4yr = sum(time_before_index>=365*4),
            n_5yr = sum(time_before_index>=365*5))

adi_window <- 30

adi_ids <- adi_visits %>% 
  filter(between(days_since_index,-adi_window,adi_window)) %>% 
  distinct(patient_id)


index_dx_dates %>% 
  inner_join(adi_ids) %>% 
  summarise(n_1yr = sum(time_before_index>=365),
            n_2yr = sum(time_before_index>=365*2),
            n_3yr = sum(time_before_index>=365*3),
            n_4yr = sum(time_before_index>=365*4),
            n_5yr = sum(time_before_index>=365*5))

adi_window <- 60

adi_ids <- adi_visits %>% 
  filter(between(days_since_index,-adi_window,adi_window)) %>% 
  distinct(patient_id)

index_dx_dates %>% 
  inner_join(adi_ids) %>% 
  summarise(n_1yr = sum(time_before_index>=365),
            n_2yr = sum(time_before_index>=365*2),
            n_3yr = sum(time_before_index>=365*3),
            n_4yr = sum(time_before_index>=365*4),
            n_5yr = sum(time_before_index>=365*5))


### Trend Plots ----------------------------------------------------------------

adi_window <- 60

obs_window <- 365*1.5

adi_ids <- adi_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  filter(between(days_since_index,-adi_window,adi_window)) %>% 
  distinct(patient_id)

include_ids <- index_dx_dates %>% 
  filter(time_before_index>=obs_window) %>% 
  distinct(patient_id)

include_ids <- include_ids  %>% 
  left_join(mutate(adi_ids,adi_patient=TRUE)) %>% 
  mutate(adi_patient = replace_na(adi_patient,FALSE))


count(include_ids,adi_patient)

ssd_counts <- ssd_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  inner_join(include_ids) %>% 
  filter(between(days_since_index,-365*2,-1)) %>% 
  group_by(adi_patient) %>% 
  count(days_since_index) %>% 
  ungroup() 

ssd_counts %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_wrap(~adi_patient,scales = "free_y")

ssd_counts %>% 
  inner_join(count(include_ids,adi_patient,name = "denom")) %>% 
  mutate(`AIDS-defining illness` = ifelse(adi_patient,"Yes","No")) %>% 
  mutate(`AIDS-defining illness` = fct_relevel(`AIDS-defining illness`,"Yes","No")) %>% 
  mutate(frac = 100*n/denom) %>% 
  mutate(period = -days_since_index) %>% 
  filter(period <= 365*1.5) %>% 
  ggplot(aes(period,frac, color = `AIDS-defining illness`)) +
  geom_line(size = .75) +
  theme_bw() +
  theme(legend.position = c(0.3, 0.7)) +
  scale_x_reverse() +
  xlab("Days Before Index HIV Diagnosis") +
  ylab("% of Patients with SSD-Related Visit") +
  labs(color = "Patients diagnosed with \n AIDS-defining illness")
ggsave("~/OneDrive - University of Iowa/grant_proposals/internal_hiv_delay/fig1.pdf",width = 4, height = 3)



### SSD Plot by condition ------------------------------------------------------

adi_window <- 60

obs_window <- 365*1.5

adi_ids <- adi_visits %>% 
  distinct(patient_id,days_since_index) %>% 
  filter(between(days_since_index,-adi_window,adi_window)) %>% 
  distinct(patient_id)

include_ids <- index_dx_dates %>% 
  filter(time_before_index>=obs_window) %>% 
  distinct(patient_id)

include_ids <- include_ids  %>% 
  left_join(mutate(adi_ids,adi_patient=TRUE)) %>% 
  mutate(adi_patient = replace_na(adi_patient,FALSE))

ssd_visits %>% 
  inner_join(include_ids) %>% 
  filter(between(days_since_index,-obs_window,-1)) %>% 
  count(adi_patient,days_since_index,label) %>% 
  inner_join(count(include_ids,adi_patient,name = "denom")) %>% 
  mutate(frac = 100*n/denom) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(`AIDS-defining illness` = ifelse(adi_patient,"Yes","No")) %>% 
  mutate(`AIDS-defining illness` = fct_relevel(`AIDS-defining illness`,"Yes","No")) %>% 
  ggplot(aes(period,frac, color = `AIDS-defining illness`)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~label, scales = "free_y") +
  theme(legend.position = c(0.8, 0.2)) +
  scale_x_reverse() +
  xlab("Days Before Index HIV Diagnosis") +
  ylab("% of Patients with SSD-Related Visit")


ssd_visits %>% 
  inner_join(include_ids) %>% 
  filter(between(days_since_index,-obs_window,-1)) %>% 
  mutate(weeks_since_index = days_since_index %/% 7) %>% 
  count(adi_patient,weeks_since_index,label) %>% 
  inner_join(count(include_ids,adi_patient,name = "denom")) %>% 
  mutate(frac = 100*n/denom) %>% 
  mutate(week = -weeks_since_index) %>% 
  mutate(`AIDS-defining illness` = ifelse(adi_patient,"Yes","No")) %>% 
  mutate(`AIDS-defining illness` = fct_relevel(`AIDS-defining illness`,"Yes","No")) %>% 
  ggplot(aes(week,frac, color = `AIDS-defining illness`)) +
  geom_line(size = 1) +
  theme_bw() +
  facet_wrap(~label, scales = "free_y") +
  theme(legend.position = "bottom") +
  scale_x_reverse() +
  xlab("Weeks Before Index HIV Diagnosis") +
  ylab("% of Patients with SSD-Related Visit") +
  labs(color = "Patients diagnosed with AIDS-defining illness")
ggsave("~/OneDrive - University of Iowa/grant_proposals/internal_hiv_delay/fig3.pdf",width = 7.5, height = 4.5)
  


### Plot with delay areas shaded -----------------------------------------------


tmp_ssd_counts <- ssd_visits %>% 
  inner_join(include_ids) %>% 
  filter(between(days_since_index,-obs_window,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  filter(days_since_index<0) %>%
  mutate(dow = days_since_index %% 7) %>%
  mutate(dow = as.factor(dow))

assemble_cp_fit(tmp_ssd_counts,cp = -200) %>% 
  filter(model == "quad") %>% 
  mutate(excess = ifelse(excess==0,NA,combined)) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_point() +
  geom_line(aes(y = excess), color = "red", size = 1.5) +
  geom_line(aes(y = expected), color = "blue", size = 1.5) +
  geom_ribbon(aes(ymin=expected,ymax = combined), fill = "red", alpha = 0.5) +
  geom_ribbon(aes(ymin=0,ymax = expected), fill = "blue", alpha = 0.5) +
  theme_bw() +
  scale_x_reverse() +
  ylim(0,400) +
  ylab("Number of Patients") +
  xlab("Days Before Index HIV Diagnosis")

trends_ssd %>%
  filter(cp ==14) %>%  
  filter(period<=75) %>% 
  select(period,"Observed"=n,"Expected"=pred,cp)  %>% 
  ggplot(aes(period,Observed)) +
  geom_point(size = 2) +
  geom_vline(aes(xintercept = cp), linetype =2) +
  # geom_vline(aes(xintercept = 1), linetype =1) +
  theme_bw() +
  geom_ribbon(aes(ymin=0,ymax = Expected-425), fill = "blue", alpha = 0.5) +
  geom_ribbon(aes(ymin=Expected,ymax = Observed), fill = "red", alpha = 0.5) +
  scale_color_manual(values = c("red","black")) +
  geom_line(aes(y = Expected), color = "blue", size = 1.2) +
  ylab("Total SSD visits") +
  xlab("Days Before Index Sepsis Diagnosis") +
  scale_x_reverse(breaks = c(1,14,75),expand = c(0.014, 0)) +
  scale_y_continuous(expand = c(0.021,0))



### Use new bootstrapping approach ---------------------------------------------

load(paste0("~/Data/Statepi_Diagnosis/prelim_results/hiv/delay_results/all_dx_visits.RData"))

upper_bound <- 365*2

boot_counts_ssd_adi <- tibble(bootstrap = 1:100) %>%
  mutate(vis_counts = map(bootstrap,
                          ~compute_boot_count(index_dates = inner_join(index_dx_dates,
                                                                       tmp_ids,
                                                                       by = join_by(patient_id)),
                                              dx_visit_data = inner_join(ssd_visits,
                                                                         tmp_ids,
                                                                         by = join_by(patient_id)))))


boot_counts_ssd_non_adi <- tibble(bootstrap = 1:100) %>%
  mutate(vis_counts = map(bootstrap,
                          ~compute_boot_count(index_dates = anti_join(index_dx_dates,
                                                                       tmp_ids,
                                                                       by = join_by(patient_id)),
                                              dx_visit_data = anti_join(ssd_visits,
                                                                         tmp_ids,
                                                                         by = join_by(patient_id)))))
tmp <- boot_counts_ssd_non_adi$vis_counts[[1]] %>% 
  distinct(days_since_index,dow)

boot_counts_ssd_adi <- boot_counts_ssd_adi %>% 
  mutate(vis_counts = map(vis_counts,~left_join(tmp,.,by = join_by(days_since_index, dow)))) %>% 
  mutate(vis_counts = map(vis_counts,~mutate(.,n = replace_na(n,0)))) 


boot_counts_ssd_non_adi <- boot_counts_ssd_non_adi %>% 
  mutate(vis_counts = map(vis_counts,~left_join(tmp,.,by = join_by(days_since_index, dow)))) %>% 
  mutate(vis_counts = map(vis_counts,~mutate(.,n = replace_na(n,0)))) 

cp_range <- -250:(-100)

boot_counts_ssd_adi <- boot_counts_ssd_adi %>%
  mutate(fits = map(vis_counts,~fit_trends(.,cp_range = cp_range))) %>%
  select(bootstrap,fits)


select(filter(ssd_counts,adi_patient),days_since_index,n)

out_of_sample_mse_adi <- boot_counts_ssd_adi %>%
  unnest(fits) %>%
  select(bootstrap,cp,model,days_since_index,combined) %>%
  left_join(select(filter(ssd_counts,adi_patient),days_since_index,n)) %>% 
  group_by(bootstrap,cp,model) %>%
  mutate(n = replace_na(n,0)) %>% 
  summarise(mse = mean((combined-n)^2)) %>%
  group_by(cp,model) %>%
  summarise(out_mse = mean(mse))

out_of_sample_mse_adi %>% 
  ungroup() %>% 
  arrange(out_mse)


