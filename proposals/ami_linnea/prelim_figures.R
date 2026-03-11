

library(tidyverse)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/ami/ami.db")

# A figure with 4 panels of SSDs (ICD9, ICD10) — cough (R05, 7862), unspecified abdominal pain (R19, 78900), 
#  shortness of breath (R0602, 78605), and chest pain (R079, 78650)
# A figure like the one that you keep using for TB, but for AMI.

### Collect Data on HPC --------------------------------------------------------

ssd_codes <- tibble(dx = c("7862","R05",
                           "78900","R19",
                           "78605","R0602",
                           "78650","R079"),
                    dx_ver = c(9L,10L,
                               9L,10L,
                               9L,10L,
                               9L,10L),
                    cond = c("Cough","Cough",
                             "Unspecified abdominal pain","Unspecified abdominal pain",
                             "Shortness of breath","Shortness of breath",
                             "Chest pain","Chest pain"))

# Collect DX visits

ssd_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(ssd_codes, copy = TRUE) %>% 
  filter(days_since_index>=-365,
         days_since_index<= -1) %>% 
  collect()


all_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(days_since_index>=-365,
         days_since_index<= -1) %>% 
  distinct(patient_id,days_since_index) %>% 
  collect()

# Collect Enrollment

collapse_enroll <- bind_rows(tbl(db,"ccae_mdcr_collapse_enroll") %>% 
                               collect(),
                             tbl(db,"medicaid_collapse_enroll") %>% 
                               collect())

# Collect Index Date

index_dates <- tbl(db,"index_dx_dates") %>% 
  collect()


### Filter DX visits 365 -------------------------------------------------------

include_ids <- collapse_enroll %>% 
  inner_join(select(index_dates,patient_id,index_date)) %>% 
  filter(dtstart<=index_date,dtend>=index_date) %>% 
  mutate(enroll_before = index_date-dtstart) %>% 
  filter(enroll_before>=365) %>% 
  distinct(patient_id)

ssd_vis <- ssd_vis %>% 
  inner_join(include_ids)

all_vis <- all_vis %>% 
  inner_join(include_ids)

# generate counts

pop_count <- nrow(include_ids)

ssd_counts <- ssd_vis %>% 
  distinct(patient_id,cond,days_since_index) %>% 
  count(cond,days_since_index)

any_ssd_counts <- ssd_vis %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index)

any_vis_counts <- all_vis %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index)


save(ssd_counts,ssd_counts,any_ssd_counts,any_vis_counts,
     file = "/Shared/AML/grant_proposal_work/ami_linnea_2026/data/vis_counts_365.RData")


### Run locally ----------------------------------------------------------------

load("/Volumes/AML/grant_proposal_work/ami_linnea_2026/data/vis_counts_365.RData")

ssd_counts %>% 
  filter(days_since_index ==0)

ssd_counts %>% 
  filter(days_since_index>= -180) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  facet_wrap(~cond, scales = "free_y") +
  theme_bw() +
  scale_x_reverse() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis")
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/4_part_ssd_180days.pdf",
       width = 5, height = 4)

ssd_counts %>% 
  filter(days_since_index>= -90) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  facet_wrap(~cond, scales = "free_y") +
  theme_bw() +
  scale_x_reverse() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis")
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/4_part_ssd_90days.pdf",
       width = 5, height = 4)



any_vis_counts %>% 
  filter(days_since_index>= -180) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis")

any_ssd_counts %>% 
  filter(days_since_index>= -180) %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis")


### Fit models -----------------------------------------------------------------

fit_model <- function(data,cp,model="lm",periodicity = FALSE, return_fit=FALSE, return_pred = FALSE){
  
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
  
  if (!return_fit & !return_pred){
    return(rmse)
  } else {
    
    out <- list(rmse = rmse)
    
    if (return_fit){
      out[["fit"]] <- fit
    }
    
    if (return_pred){
      
      pred <- data %>% 
        mutate(before_cp = TRUE) %>%
        mutate(pred1=predict(fit,newdata=.)) %>%
        mutate(before_cp = period>=cp) %>% 
        mutate(pred2=predict(fit,newdata=.)) %>% 
        select(period,n,pred1,pred2)
      
      out[["pred"]] <- pred
    }
    
    return(out)
    
  }
}


tmp <- any_ssd_counts %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  filter(period<=180)


fit_data <- fit_model(data = tmp, model = "quad",cp = 14,return_pred = TRUE)

fit_data$pred %>% 
  mutate(obs = ifelse(period>14, NA, n)) %>% 
  ggplot(aes(period,n)) +
  scale_x_reverse() +
  geom_line(size = 1) +
  geom_line(aes(period,pred1), color = "red") +
  geom_ribbon(aes(ymin = 0,ymax = pred1), fill = "red", alpha = 0.25) +
  geom_ribbon(aes(ymin = pred1,ymax = obs), fill = "blue", alpha = 0.25) +
  theme_bw() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis") +
  geom_vline(aes(xintercept = 14), linetype = 2)
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/excess_visits_14cp_180days.pdf",
       width = 4, height = 3)

fit_data$pred %>% 
  filter(period<=90) %>% 
  mutate(obs = ifelse(period>14, NA, n)) %>% 
  ggplot(aes(period,n)) +
  scale_x_reverse() +
  geom_line(size = 1) +
  geom_line(aes(period,pred1), color = "red") +
  geom_ribbon(aes(ymin = 0,ymax = pred1), fill = "red", alpha = 0.25) +
  geom_ribbon(aes(ymin = pred1,ymax = obs), fill = "blue", alpha = 0.25) +
  theme_bw() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis") +
  geom_vline(aes(xintercept = 14), linetype = 2)
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/excess_visits_14cp_90days.pdf",
       width = 4, height = 3)

fit_data <- fit_model(data = tmp, model = "quad",cp = 21,return_pred = TRUE)

fit_data$pred %>% 
  mutate(obs = ifelse(period>21, NA, n)) %>% 
  ggplot(aes(period,n)) +
  scale_x_reverse() +
  geom_line(size = 1) +
  geom_line(aes(period,pred1), color = "red") +
  geom_ribbon(aes(ymin = 0,ymax = pred1), fill = "red", alpha = 0.25) +
  geom_ribbon(aes(ymin = pred1,ymax = obs), fill = "blue", alpha = 0.25) +
  theme_bw() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis") +
  geom_vline(aes(xintercept = 21), linetype = 2)
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/excess_visits_21cp_180days.pdf",
       width = 4, height = 3)

fit_data$pred %>% 
  filter(period<=90) %>% 
  mutate(obs = ifelse(period>21, NA, n)) %>% 
  ggplot(aes(period,n)) +
  scale_x_reverse() +
  geom_line(size = 1) +
  geom_line(aes(period,pred1), color = "red") +
  geom_ribbon(aes(ymin = 0,ymax = pred1), fill = "red", alpha = 0.25) +
  geom_ribbon(aes(ymin = pred1,ymax = obs), fill = "blue", alpha = 0.25) +
  theme_bw() +
  ylab("Number of patient visits") +
  xlab("Days before AMI diagnosis") +
  geom_vline(aes(xintercept = 21), linetype = 2)
ggsave("~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/figures/excess_visits_21cp_90days.pdf",
       width = 4, height = 3)


#############
#### New ####



library(tidyverse)
library(codeBuildr)

db <- src_sqlite("/Shared/AML/truven_extracts/small_dbs/ami/ami.db")

ssd_codes <- bind_rows(bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe("280")),
                    tibble(dx_ver = 10L,
                           dx = children_safe("D50"))) %>% 
            mutate(name = "Iron Deficiency Anemia",
                   group = "unrelated"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe("585")),
                    tibble(dx_ver = 10L,
                           dx = children_safe("N18"))) %>% 
            mutate(name = "Chronic Kidney Disease",
                   group = "unrelated"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = c(children_safe("690"),children_safe("691"),children_safe("692"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("L230",
                                                "L231",
                                                "L232",
                                                "L233",
                                                "L234",
                                                "L235",
                                                "L236",
                                                "L237",
                                                "L238",
                                                "L239")) )) %>% 
            mutate(name = "Atopic Dermatitis",
                   group = "unrelated"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe(c("250"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("E10")) )) %>% 
            mutate(name = "Type 1 Diabetes",
                   group = "unrelated"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe(c("486"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("J189")) )) %>% 
            mutate(name = "Pneumonia",
                   group = "ssds"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe(c("4660"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("J209")) )) %>% 
            mutate(name = "Acute Bronchitis",
                   group = "ssds"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe(c("49390"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("J459")) )) %>% 
            mutate(name = "Asthma",
                   group = "ssds"),
          
          bind_rows(tibble(dx_ver = 9L,
                           dx = children_safe(c("53081"))),
                    tibble(dx_ver = 10L,
                           dx = children_safe(c("K219")) )) %>% 
            mutate(name = "Reflux",
                   group = "ssds")
          
          ) 
  


ssd_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(ssd_codes, copy = TRUE) %>% 
  filter(days_since_index>=-365,
         days_since_index<= -1) %>% 
  collect()

# Fix Diabetes

dx_visits <- db %>% 
  tbl("all_dx_visits") %>% 
  filter(days_since_index>=-365,
         days_since_index<= -1) %>% 
  collect()

tmp1 <- dx_visits %>% 
  filter(str_detect(dx,"^E10"),
         dx_ver==10) %>% 
  mutate(name = "Type 1 Diabetes",
         group = "unrelated")

tmp2 <- dx_visits %>% 
  filter(str_detect(dx, "^250[0-9][13]"),
         dx_ver==9) %>% 
  mutate(name = "Type 1 Diabetes",
         group = "unrelated")

diabetes_vis <- bind_rows(tmp1,tmp2)

ssd_vis <- bind_rows(ssd_vis %>% 
                       filter(name!="Type 1 Diabetes"),
                     diabetes_vis)


index_dates <- tbl(db,"index_dx_dates") %>% 
  collect()


ssd_counts <- ssd_vis %>% 
  distinct(patient_id,name,group,days_since_index) %>% 
  count(name,group,days_since_index)

ssd_counts2 <- ssd_counts

save(ssd_counts2,
     file = "/Shared/AML/grant_proposal_work/ami_linnea_2026/data/ssd_unrelated_counts_365.RData")

rm(list = ls())

## Locally 

load("/Volumes/AML/grant_proposal_work/ami_linnea_2026/data/ssd_unrelated_counts_365.RData")


load("/Volumes/AML/grant_proposal_work/ami_linnea_2026/data/vis_counts_365.RData")

plot_data <- bind_rows(mutate(ssd_counts,
                 group = "set1") %>% 
            rename(name = cond),
          ssd_counts2)

write_csv(plot_data,
          "~/OneDrive - University of Iowa/grant_proposals/ami_linnea_2026/new_plots/plot_data.csv")

rm(list = ls())


