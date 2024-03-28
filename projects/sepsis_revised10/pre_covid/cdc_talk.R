
library(tidyverse)

load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/index_cases.RData")

index_cases

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/sepsis_revised10/sepsis_revised10.db")

dx_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(between(days_since_index,-180,-1)) %>% 
  filter(dx_ver==10) %>% 
  collect()

dx_vis <- index_cases %>% 
  select(patient_id) %>% 
  inner_join(dx_vis)

dx_vis %>% 
  filter(dx=="R0602") %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(xintercept = 14, linetype = 2) +
  ylab("Number of Visits") +
  xlab("Days Before Index Diagnosis") +
  ggtitle("R0602 - Shortness of Breath")

dx_vis %>% 
  filter(dx=="R509") %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(xintercept = 14, linetype = 2) +
  ylab("Number of Visits") +
  xlab("Days Before Index Diagnosis") +
  ggtitle("R509 - Fever, unspecified")


# R4182 − Altered mental status, unspecified

dx_vis %>% 
  filter(dx=="R4182") %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(xintercept = 14, linetype = 2) +
  ylab("Number of Visits") +
  xlab("Days Before Index Diagnosis") +
  ggtitle("R4182 - Altered mental status, unspecified")


# R52 − Pain, unspecified

dx_vis %>% 
  filter(dx=="R52") %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(xintercept = 14, linetype = 2) +
  ylab("Number of Visits") +
  xlab("Days Before Index Diagnosis") +
  ggtitle("R52 - Pain, unspecified")

#### Compute Trends ####

ssd_codes <- codeBuildr::load_ssd_codes("sepsis_revised10") %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver) %>%
  distinct() %>%
  filter(!is.na(dx))

ssd_count <- ssd_codes %>% 
  select(dx) %>% 
  inner_join(dx_vis) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index)

ssd_count %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(xintercept = 14, linetype = 2) +
  ylab("Number of Visits") +
  xlab("Days Before Index Diagnosis") +
  ggtitle("R52 - Pain, unspecified")

ssd_vis_count <- ssd_count %>% 
  mutate(period = - days_since_index) %>% 
  filter(period>0) %>% 
  mutate(dow = as.factor(period %% 7))

i <- 7


fit1_ssd <- lm(n ~ period + dow, filter(ssd_vis_count, period>i))
fit2_ssd <- lm(n ~ log(period) + dow, filter(ssd_vis_count, period>i))
fit3_ssd <- lm(n ~ poly(period,degree = 2) + dow, filter(ssd_vis_count, period>i))
fit4_ssd <- lm(n ~ poly(period,degree = 3) + dow, filter(ssd_vis_count, period>i))

ssd_vis_count %>% 
  mutate(Linear = predict(fit1_ssd,newdata =.)) %>% 
  mutate(Exponential = predict(fit2_ssd,newdata =.)) %>% 
  mutate(Quadratic = predict(fit3_ssd,newdata =.)) %>% 
  mutate(Cubic = predict(fit4_ssd,newdata =.)) %>% 
  select(period,n,Linear:Cubic) %>% 
  gather(key = Model, value = value, -period, -n) %>% 
  ggplot(aes(period,value, color = Model)) +
  geom_line() +
  geom_line(aes(period,n), color = "black") + 
  geom_vline(xintercept = 14, linetype = 2) +
  scale_x_reverse() +
  theme_minimal() +
  scale_color_manual(values = c("orange","red","purple","blue")) +
  xlab("Days Before Sepsis Diagnosis") +
  ylab("Number of Visits for Signs or Symptoms")


ssd_vis_count %>% 
  ggplot(aes(period,n)) +
  geom_line() +
  geom_vline(xintercept = 14, linetype = 2) +
  scale_x_reverse() +
  theme_minimal() +
  xlab("Days Before Sepsis Diagnosis") +
  ylab("Number of Visits for Signs or Symptoms")


rm(list = ls())
