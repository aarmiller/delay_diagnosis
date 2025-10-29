

library(tidyverse)

load("/Volumes/Statepi_Diagnosis/prelim_results/cocci/delay_results/all_dx_visits.RData")


condition <- "cocci"

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

ssd_fits <- tibble(model = c("lm","quad","cubic","exp")) %>% 
  mutate(label = c("Linear", "Quadratic", "Cubic", "Exponential")) %>% 
  mutate(periodicity = map(model,~c(T,F))) %>% 
  unnest(periodicity) %>% 
  mutate(label = ifelse(periodicity==TRUE,paste0(label," w/ periodicity"),label)) %>% 
  mutate(cp_res = map2(model,periodicity,
                       ~fit_cp_range(count_data_ssd,
                                     cp_range = seq(from = 7, to = 300, by = 7),
                                     model = .x,
                                     periodicity = .y)))

ssd_fits <- ssd_fits %>% 
  unnest(cp_res) %>% 
  mutate(rmse = map_dbl(fit_res,~.$rmse)) %>% 
  mutate(pred = map(fit_res,~.$pred)) %>% 
  select(-fit_res)

ssd_fits <- ssd_fits %>% 
  mutate(miss = map2(pred,cp,~mutate(.x,cp = .y))) %>% 
  mutate(miss = map(miss,~mutate(., n_miss = ifelse(n>pred1,n-pred1,0)))) %>% 
  mutate(miss = map(miss, ~filter(., period<=cp))) %>% 
  mutate(miss = map(miss,~summarise(.,n_miss = sum(n_miss)))) %>% 
  unnest(miss)

appendix_figure1 <- ssd_fits %>% 
  filter(model == "quad",
         periodicity == TRUE) %>% 
  filter(cp %in% c(63,70,77,
                   84,91,98,
                   105,112,119)) %>% 
  mutate(weeks = cp/7) %>% 
  mutate(label = paste0(weeks,"-Week Opportunity Window")) %>% 
  unnest(pred) %>% 
  filter(weeks %in% 12:14) %>%
  ggplot(aes(period,n)) +
  geom_line(aes(y = pred2), color = "blue", size = .75) +
  geom_line(aes(y = pred1), color ="red", size = .75) +
  geom_point(alpha = 0.5) +
  scale_x_reverse() +
  facet_wrap(~label) + 
  geom_vline(aes(xintercept = cp), linetype = 2) +
  theme_bw() +
  ylab("Number of SSD Visits") +
  xlab("Days Before Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/submissions/OFID/revisions/appendix_figure1.pdf",
       width = 10, height = 5,dpi = 600,units = "in",
       plot = appendix_figure1 )
