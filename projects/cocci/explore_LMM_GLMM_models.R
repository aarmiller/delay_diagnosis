

#### Delay Patient All ---------------------------------------------------------
library(lme4)

get_delay_res_LMM <- function(trial_val){
  
  
  tmp1 <- full_reg_data_dur %>% filter(trial==trial_val)
  
  reg_data <- tmp1 %>% 
    select(boot_sample) %>% 
    unnest(boot_sample) %>% 
    select(patient_id) %>% 
    inner_join(reg_demo, by = "patient_id") %>% 
    inner_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # left_join(tmp1 %>% select(data) %>% unnest(data), by = "patient_id") %>%
    # mutate(duration = replace_na(duration,0L)) %>% 
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(duration, age_cat, female, rural, source, year, month,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp,
           #top2_high_inc_state_baddley,
           state_abb,
           resp_antibiotic_drugs_window, inhalers_window) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  fit <- lme4::lmer(duration ~ age_cat + female + rural + source + year + month + 
               asthma_prior_cp + copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp + 
               resp_antibiotic_drugs_window + inhalers_window +(1|state_abb),
             data = reg_data)
  
  RE_temp <- lme4::ranef(fit)$state_abb
  FE_temp <- lme4::fixef(fit)
  
  RE <- tibble(intercept = RE_temp[,1],
               state_abb = rownames(RE_temp))
  FE <- tibble(term = names(FE_temp),
               estimate = FE_temp)
  
  return(list(RE = RE,
              FE = FE))
         
}


miss_delay_pat_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                         function(x){get_delay_res_LMM(x)}, 
                                         mc.cores = num_cores)

RE <- do.call("rbind", lapply(miss_delay_pat_res, function(x){x$RE} ))
FE <- do.call("rbind", lapply(miss_delay_pat_res, function(x){x$FE} ))


miss_delay_dur_LMM_FE <- FE %>% 
  group_by(term) %>% 
    summarise(est = median(estimate, na.rm = T),
              low = quantile(estimate, probs = c(0.025), na.rm = T),
              high = quantile(estimate, probs = c(0.975), na.rm = T),
              low_90 = quantile(estimate, probs = c(0.05), na.rm = T),
              high_90 = quantile(estimate, probs = c(0.95), na.rm = T) )

miss_delay_dur_LMM_RE <- RE %>% 
  group_by(state_abb) %>% 
  summarise(intercept_est = median(intercept, na.rm = T),
            low = quantile(intercept, probs = c(0.025), na.rm = T),
            high = quantile(intercept, probs = c(0.975), na.rm = T),
            low_90 = quantile(intercept, probs = c(0.05), na.rm = T),
            high_90 = quantile(intercept, probs = c(0.95), na.rm = T))

save(miss_delay_dur_LMM_FE, miss_delay_dur_LMM_RE, file = paste0(out_path,"LMM_state_data.RData"))
print("LMM finished")

#### Miss opp All ---------------------------------------------------------

get_miss_opp_res_GLMM <- function(trial_val){
  
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  reg_data <- bind_rows(tmp1 %>% select(data) %>% 
                          unnest(data) %>% 
                          mutate(miss=TRUE),
                        tmp1 %>% select(boot_sample) %>% 
                          unnest(boot_sample) %>% 
                          mutate(miss=FALSE)) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% #remove year and month as that is for index date in reg_demo
    filter(year > 2001) %>% 
    mutate(year = as.factor(year),
           month = factor(month, levels = 1:12)) %>% 
    select(miss, inpatient, age_cat, female, rural, source, weekend, year, month,
           asthma_prior_cp, copd_prior_cp, chest_ct_prior_cp, chest_xray_prior_cp, 
           #top2_high_inc_state_baddley
           state_abb) %>% 
    drop_na() # have to do at the end as obs not in bootsample and boot_id not in data
  
  
  fit <- lme4::glmer(miss ~ inpatient + age_cat + female + rural + source + weekend + year + month + 
               asthma_prior_cp + copd_prior_cp + chest_ct_prior_cp + chest_xray_prior_cp +(1|state_abb),
             data = reg_data,
             family = binomial(link = "logit"))
  
  RE_temp <- lme4::ranef(fit)$state_abb
  FE_temp <- lme4::fixef(fit)
  
  RE <- tibble(intercept = RE_temp[,1],
               state_abb = rownames(RE_temp))
  FE <- tibble(term = names(FE_temp),
               estimate = FE_temp)
  
  return(list(RE = RE,
              FE = FE))
  
}


miss_opp_res <- parallel::mclapply(1:max(full_reg_data$trial),
                                         function(x){get_miss_opp_res_GLMM(x)}, 
                                         mc.cores = num_cores)

RE <- do.call("rbind", lapply(miss_opp_res, function(x){x$RE} ))
FE <- do.call("rbind", lapply(miss_opp_res, function(x){x$FE} ))


miss_opp_GLMM_FE <- FE %>% 
  group_by(term) %>% 
  summarise(est = median(exp(estimate), na.rm = T),
            low = quantile(exp(estimate), probs = c(0.025), na.rm = T),
            high = quantile(exp(estimate), probs = c(0.975), na.rm = T),
            low_90 = quantile(exp(estimate), probs = c(0.05), na.rm = T),
            high_90 = quantile(exp(estimate), probs = c(0.95), na.rm = T) )

miss_opp_GLMM_RE <- RE %>% 
  group_by(state_abb) %>% 
  summarise(intercept_est = median(intercept, na.rm = T),
            low = quantile(intercept, probs = c(0.025), na.rm = T),
            high = quantile(intercept, probs = c(0.975), na.rm = T),
            low_90 = quantile(intercept, probs = c(0.05), na.rm = T),
            high_90 = quantile(intercept, probs = c(0.95), na.rm = T))

save(miss_opp_GLMM_FE, miss_opp_GLMM_RE, file = paste0(out_path,"GLMM_state_data.RData"))
print("GLMM finished")

################################################################################
# 
# library(tidyverse)
# library(usmap)
# out_path <- "/Volumes/Statepi_Diagnosis/projects/cocci/risk_models/"
# load(paste0(out_path,"LMM_state_data.RData"))
# 
# states_w_info_LMM <- usmap::statepop %>% select(-pop_2015) %>%
#   inner_join(miss_delay_dur_LMM_RE %>% rename(abbr = state_abb))
# 
# plot_usmap(data = states_w_info_LMM, values = "intercept_est")+
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", name =  "",
#                         label = scales::comma,) +
#   labs(title = "LMM for Duration of Delay", subtitle = "Predicted State-level Random Intercept") +
#   theme(legend.position = "right")
# 
# 
# 
# load(paste0(out_path,"GLMM_state_data.RData"))
# 
# states_w_info_GLMM <- usmap::statepop %>% select(-pop_2015) %>%
#   inner_join(miss_opp_GLMM_RE %>% rename(abbr = state_abb))
# 
# plot_usmap(data = states_w_info_GLMM, values = "intercept_est")+
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", name =  "",
#                        label = scales::comma,) +
#   labs(title = "GLMM for Log-odds of Miss Opportunities", subtitle = "Predicted State-level Random Intercept") +
#   theme(legend.position = "right")
# 
# 
