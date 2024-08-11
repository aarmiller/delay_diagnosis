ssd_state_fun <- function(trial_val){
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    mutate(miss=TRUE) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    group_by(location, state_name, state_abb) %>% 
    summarise(n = n()) %>%
    ungroup() %>% 
    mutate(percent = n/sum(n)*100) %>% 
    arrange(desc(n))
   
  reg_demo %>% distinct(location,state_name,state_abb ) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(n, percent),~replace_na(.,0L))
}

ssd_state <- parallel::mclapply(1:max(full_reg_data$trial),
                                           function(x){ssd_state_fun(x)}, 
                                           mc.cores = num_cores)


ssd_state <- bind_rows(ssd_state) %>% 
  group_by(location, state_name, state_abb) %>% 
  summarise(n_mean = mean(n, na.rm = T),
            n_low = quantile(n, probs = c(0.025), na.rm = T),
            n_high = quantile(n, probs = c(0.975), na.rm = T),
            percent_mean = mean(percent, na.rm = T),
            percent_low = quantile(percent, probs = c(0.025), na.rm = T),
            percent_high = quantile(percent, probs = c(0.975), na.rm = T))

ssd_state <- ssd_state %>% 
  ungroup() %>% 
  arrange(desc(n_mean)) %>% 
  mutate(nice_n = paste0(format(round(n_mean,2), nsmall =2),
                         " (",format(round(n_low,2), nsmall =2), ", ",
                         format(round(n_high,2), nsmall =2), ")"),
         nice_percent = paste0(format(round(percent_mean,2), nsmall =2),
                               " (",format(round(percent_low,2), nsmall =2), ", ",
                               format(round(percent_high,2), nsmall =2), ")")) 

write_csv(ssd_state, file = paste0(out_path,"ssd_counts_by_state.csv"))  

################################################################################
ssd_baddley_fun <- function(trial_val){
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    mutate(miss=TRUE) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    count(top2_high_inc_state_baddley) %>% 
    mutate(percent = n/sum(n)*100) %>% 
    arrange(desc(n))
  
  reg_demo %>% distinct(top2_high_inc_state_baddley) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(n, percent),~replace_na(.,0L))
}

ssd_baddley <- parallel::mclapply(1:max(full_reg_data$trial),
                                function(x){ssd_baddley_fun(x)}, 
                                mc.cores = num_cores)


ssd_baddley <- bind_rows(ssd_baddley) %>% 
  group_by(top2_high_inc_state_baddley) %>% 
  summarise(n_mean = mean(n, na.rm = T),
            n_low = quantile(n, probs = c(0.025), na.rm = T),
            n_high = quantile(n, probs = c(0.975), na.rm = T),
            percent_mean = mean(percent, na.rm = T),
            percent_low = quantile(percent, probs = c(0.025), na.rm = T),
            percent_high = quantile(percent, probs = c(0.975), na.rm = T))

ssd_baddley <- ssd_baddley %>% 
  ungroup() %>% 
  arrange(desc(n_mean)) %>% 
  mutate(nice_n = paste0(format(round(n_mean,2), nsmall =2),
                         " (",format(round(n_low,2), nsmall =2), ", ",
                         format(round(n_high,2), nsmall =2), ")"),
         nice_percent = paste0(format(round(percent_mean,2), nsmall =2),
                               " (",format(round(percent_low,2), nsmall =2), ", ",
                               format(round(percent_high,2), nsmall =2), ")")) 

write_csv(ssd_baddley, file = paste0(out_path,"ssd_counts_baddley.csv"))  


################################################################################

ssd_west_state_fun <- function(trial_val){
  
  west_states <- c("WA", "OR", "CA",
                   "ID", "NV", "AZ",
                   "MT", "WY", "CO",
                   "UT", "NM")
  
  tmp1 <- full_reg_data %>% 
    filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    mutate(miss=TRUE) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    mutate(west_states = ifelse(is.na(state_abb), NA, (state_abb %in% west_states)*1.0)) %>%
    count(west_states) %>% 
    mutate(percent = n/sum(n)*100) %>% 
    arrange(desc(n))
  
   tibble(west_states = c(0,1,NA)) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(n, percent),~replace_na(.,0L))
}

ssd_west_state <- parallel::mclapply(1:max(full_reg_data$trial),
                                  function(x){ssd_west_state_fun(x)}, 
                                  mc.cores = num_cores)


ssd_west_state <- bind_rows(ssd_west_state) %>% 
  group_by(west_states) %>% 
  summarise(n_mean = mean(n, na.rm = T),
            n_low = quantile(n, probs = c(0.025), na.rm = T),
            n_high = quantile(n, probs = c(0.975), na.rm = T),
            percent_mean = mean(percent, na.rm = T),
            percent_low = quantile(percent, probs = c(0.025), na.rm = T),
            percent_high = quantile(percent, probs = c(0.975), na.rm = T))

ssd_west_state <- ssd_west_state %>% 
  ungroup() %>% 
  arrange(desc(n_mean)) %>% 
  mutate(nice_n = paste0(format(round(n_mean,2), nsmall =2),
                         " (",format(round(n_low,2), nsmall =2), ", ",
                         format(round(n_high,2), nsmall =2), ")"),
         nice_percent = paste0(format(round(percent_mean,2), nsmall =2),
                               " (",format(round(percent_low,2), nsmall =2), ", ",
                               format(round(percent_high,2), nsmall =2), ")"))

write_csv(ssd_west_state, file = paste0(out_path,"ssd_counts_west_states.csv"))  


################################################################################
####################
## Delay Duration ##
####################

full_reg_data_dur1 <- full_reg_data %>% 
  mutate(data = map(data, ~inner_join(., sim_obs_reduced, by = c("obs", "patient_id")))) 

ssd_state_fun <- function(trial_val){
  
  tmp1 <- full_reg_data_dur1 %>% filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    group_by(location, state_name, state_abb) %>% 
    summarise(duration_max = -min(days_since_index),
              duration_min = -max(days_since_index),
              duration_median = median(-days_since_index),
              duration_q1 = quantile(-days_since_index, probs = 0.25),
              duration_q3 = quantile(-days_since_index, probs = 0.75)) %>% 
    ungroup() %>% 
    arrange(desc(duration_median))
  
  reg_demo %>% distinct(location,state_name,state_abb ) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(duration_max, duration_min, duration_median, duration_q1, duration_q3),~replace_na(.,0L))
  
}

ssd_state <- parallel::mclapply(1:max(full_reg_data$trial),
                                function(x){ssd_state_fun(x)}, 
                                mc.cores = num_cores)


ssd_state <- bind_rows(ssd_state) %>% 
  group_by(location, state_name, state_abb) %>% 
  summarise(ssd_time_max = mean(duration_max, na.rm = T),
            ssd_time_max_low = quantile(duration_max, probs = c(0.025), na.rm = T),
            ssd_time_max_high = quantile(duration_max, probs = c(0.975), na.rm = T),
            ssd_time_min = mean(duration_min, na.rm = T),
            ssd_time_min_low = quantile(duration_min, probs = c(0.025), na.rm = T),
            ssd_time_min_high = quantile(duration_min, probs = c(0.975), na.rm = T),
            ssd_time_q1 = mean(duration_q1, na.rm = T),
            ssd_time_q1_low = quantile(duration_q1, probs = c(0.025), na.rm = T),
            ssd_time_q1_high = quantile(duration_q1, probs = c(0.975), na.rm = T),
            ssd_time_q3 = mean(duration_q3, na.rm = T),
            ssd_time_q3_low = quantile(duration_q3, probs = c(0.025), na.rm = T),
            ssd_time_q3_high = quantile(duration_q3, probs = c(0.975), na.rm = T),
            ssd_time_median = mean(duration_median, na.rm = T),
            ssd_time_median_low = quantile(duration_median, probs = c(0.025), na.rm = T),
            ssd_time_median_high = quantile(duration_median, probs = c(0.975), na.rm = T)) %>% 
  arrange(desc(ssd_time_median))

write_csv(ssd_state, file = paste0(out_path,"ssd_duration_by_state.csv"))  


################################################################################

ssd_baddley_fun <- function(trial_val){
  
  tmp1 <- full_reg_data_dur1 %>% filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    group_by(top2_high_inc_state_baddley) %>% 
    summarise(duration_max = -min(days_since_index),
              duration_min = -max(days_since_index),
              duration_median = median(-days_since_index),
              duration_q1 = quantile(-days_since_index, probs = 0.25),
              duration_q3 = quantile(-days_since_index, probs = 0.75)) %>% 
    ungroup() %>% 
    arrange(desc(duration_median))
  
  reg_demo %>% distinct(top2_high_inc_state_baddley) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(duration_max, duration_min, duration_median, duration_q1, duration_q3),~replace_na(.,0L))
}

ssd_baddley <- parallel::mclapply(1:max(full_reg_data$trial),
                                  function(x){ssd_baddley_fun(x)}, 
                                  mc.cores = num_cores)


ssd_baddley <- bind_rows(ssd_baddley) %>% 
  group_by(top2_high_inc_state_baddley) %>% 
  summarise(ssd_time_max = mean(duration_max, na.rm = T),
            ssd_time_max_low = quantile(duration_max, probs = c(0.025), na.rm = T),
            ssd_time_max_high = quantile(duration_max, probs = c(0.975), na.rm = T),
            ssd_time_min = mean(duration_min, na.rm = T),
            ssd_time_min_low = quantile(duration_min, probs = c(0.025), na.rm = T),
            ssd_time_min_high = quantile(duration_min, probs = c(0.975), na.rm = T),
            ssd_time_q1 = mean(duration_q1, na.rm = T),
            ssd_time_q1_low = quantile(duration_q1, probs = c(0.025), na.rm = T),
            ssd_time_q1_high = quantile(duration_q1, probs = c(0.975), na.rm = T),
            ssd_time_q3 = mean(duration_q3, na.rm = T),
            ssd_time_q3_low = quantile(duration_q3, probs = c(0.025), na.rm = T),
            ssd_time_q3_high = quantile(duration_q3, probs = c(0.975), na.rm = T),
            ssd_time_median = mean(duration_median, na.rm = T),
            ssd_time_median_low = quantile(duration_median, probs = c(0.025), na.rm = T),
            ssd_time_median_high = quantile(duration_median, probs = c(0.975), na.rm = T)) %>% 
  arrange(desc(ssd_time_median))

write_csv(ssd_baddley, file = paste0(out_path,"ssd_duration_baddley.csv"))  



################################################################################

ssd_west_fun <- function(trial_val){
  
  west_states <- c("WA", "OR", "CA",
                   "ID", "NV", "AZ",
                   "MT", "WY", "CO",
                   "UT", "NM")
  
  tmp1 <- full_reg_data_dur1 %>% filter(trial==trial_val)
  
  temp2 <- tmp1 %>% select(data) %>% 
    unnest(data) %>% 
    inner_join(reg_demo %>% select(-year, -month), by = "patient_id") %>% 
    mutate(west_states = ifelse(is.na(state_abb), NA, (state_abb %in% west_states)*1.0)) %>%
    group_by(west_states) %>% 
    summarise(duration_max = -min(days_since_index),
              duration_min = -max(days_since_index),
              duration_median = median(-days_since_index),
              duration_q1 = quantile(-days_since_index, probs = 0.25),
              duration_q3 = quantile(-days_since_index, probs = 0.75)) %>% 
    ungroup() %>% 
    arrange(desc(duration_median))
  
  tibble(west_states = c(0,1,NA)) %>% 
    left_join(temp2) %>% 
    mutate_at(vars(duration_max, duration_min, duration_median, duration_q1, duration_q3),~replace_na(.,0L))
}

ssd_west <- parallel::mclapply(1:max(full_reg_data$trial),
                                function(x){ssd_west_fun(x)}, 
                                mc.cores = num_cores)


ssd_west <- bind_rows(ssd_west) %>% 
  group_by(west_states) %>% 
  summarise(ssd_time_max = mean(duration_max, na.rm = T),
            ssd_time_max_low = quantile(duration_max, probs = c(0.025), na.rm = T),
            ssd_time_max_high = quantile(duration_max, probs = c(0.975), na.rm = T),
            ssd_time_min = mean(duration_min, na.rm = T),
            ssd_time_min_low = quantile(duration_min, probs = c(0.025), na.rm = T),
            ssd_time_min_high = quantile(duration_min, probs = c(0.975), na.rm = T),
            ssd_time_q1 = mean(duration_q1, na.rm = T),
            ssd_time_q1_low = quantile(duration_q1, probs = c(0.025), na.rm = T),
            ssd_time_q1_high = quantile(duration_q1, probs = c(0.975), na.rm = T),
            ssd_time_q3 = mean(duration_q3, na.rm = T),
            ssd_time_q3_low = quantile(duration_q3, probs = c(0.025), na.rm = T),
            ssd_time_q3_high = quantile(duration_q3, probs = c(0.975), na.rm = T),
            ssd_time_median = mean(duration_median, na.rm = T),
            ssd_time_median_low = quantile(duration_median, probs = c(0.025), na.rm = T),
            ssd_time_median_high = quantile(duration_median, probs = c(0.975), na.rm = T)) %>% 
  arrange(desc(ssd_time_median))

write_csv(ssd_west, file = paste0(out_path,"ssd_duration_west_states.csv"))  



################################################################################
# library(tidyverse)
# 
# out_path <- "/Volumes/Statepi_Diagnosis/projects/cocci/risk_models/"
# 
# # State
# ssd_count_state <- read_csv(file = paste0(out_path,"ssd_counts_by_state.csv"))  
# index_count_state <- read_csv(file = paste0(out_path,"index_counts_by_state.csv"))  
# 
# data <- index_count_state %>% mutate(group = "Index") %>%
#   bind_rows(ssd_count_state %>% mutate(group = "SSD") %>% 
#               rename(percent = percent_mean))
# 
# # wont sum to 100% by group as excluding NA
# data %>% filter(!is.na(state_abb)) %>% 
#   ggplot(aes(x = state_abb, y = percent, fill = group))+
#     geom_bar(stat="identity", position=position_dodge())+
#   ylab("Percent")+
#   xlab("State")+
#   ggtitle("Percent of visits by state")
# 
# # Baddley
# ssd_count_state <- read_csv(file = paste0(out_path,"ssd_counts_baddley.csv"))  
# index_count_state <- read_csv(file = paste0(out_path,"index_top2_high_inc_state_baddley_count.csv"))  
# 
# data <- index_count_state %>% mutate(group = "Index") %>%
#   bind_rows(ssd_count_state %>% mutate(group = "SSD") %>% 
#               rename(percent = percent_mean))
# 
# # wont sum to 100% by group as excluding NA
# data %>% filter(!is.na(top2_high_inc_state_baddley)) %>% 
#   mutate(new_lab=ifelse(top2_high_inc_state_baddley == 1, "yes", "no")) %>% 
#   ggplot(aes(x = new_lab, y = percent, fill = group))+
#   geom_bar(stat="identity", position=position_dodge())+
#   ylab("Percent")+
#   xlab("Top 2 highest incidence Baddley et al.")+
#   ggtitle("Percent of visits by Baddley et al. def")
# 
# 
# # west
# ssd_count_state <- read_csv(file = paste0(out_path,"ssd_counts_west_states.csv"))  
# index_count_state <- read_csv(file = paste0(out_path,"index_west_states_count.csv"))  
# 
# data <- index_count_state %>% mutate(group = "Index") %>%
#   bind_rows(ssd_count_state %>% mutate(group = "SSD") %>% 
#               rename(percent = percent_mean))
# 
# # wont sum to 100% by group as excluding NA
# data %>% filter(!is.na(west_states)) %>% 
#   mutate(new_lab=ifelse(west_states == 1, "yes", "no")) %>% 
#   ggplot(aes(x = new_lab, y = percent, fill = group))+
#   geom_bar(stat="identity", position=position_dodge())+
#   ylab("Percent")+
#   xlab("West States")+
#   ggtitle("Percent of visits by in west states vs. all others")
# 
# ################################################################################
# 
# # SSD time state
# ssd_time <- read_csv(file = paste0(out_path,"ssd_duration_by_state.csv"))  
# 
# data <- ssd_time %>% 
#   filter(!is.na(state_abb)) %>% 
#   select(state_abb, ssd_time_min, ssd_time_q1, ssd_time_median, ssd_time_q3, ssd_time_max) %>% 
#   pivot_longer(cols = ssd_time_min:ssd_time_max, values_to = "y", names_to = "group") %>% 
#   mutate(group = str_remove_all(group, "ssd_time_"))
# 
# data %>% 
#   mutate(group = factor(group, levels = c("min", "q1", "median", "q3", "max"))) %>% 
#   ggplot(aes(x = state_abb, y = y, fill = group))+
#   geom_bar(stat="identity", position=position_dodge()) + 
#   facet_grid(group~., scales = "free") +
#   ylab("SSD time before index")+
#   xlab("State")+
#   ggtitle("SSD time before index by state")+ 
#   guides(fill="none")
# 
# # SSD time baddley
# ssd_time <- read_csv(file = paste0(out_path,"ssd_duration_baddley.csv"))  
# 
# data <- ssd_time %>%
#   filter(!is.na(top2_high_inc_state_baddley)) %>% 
#   select(top2_high_inc_state_baddley, ssd_time_min, ssd_time_q1, ssd_time_median, ssd_time_q3, ssd_time_max) %>% 
#   pivot_longer(cols = ssd_time_min:ssd_time_max, values_to = "y", names_to = "group") %>% 
#   mutate(group = str_remove_all(group, "ssd_time_")) %>% 
#   mutate(new_lab=ifelse(top2_high_inc_state_baddley == 1, "yes", "no"))
# 
# data %>% 
#   mutate(group = factor(group, levels = c("min", "q1", "median", "q3", "max"))) %>% 
#   ggplot(aes(x = new_lab, y = y, fill = group))+
#   geom_bar(stat="identity", position=position_dodge()) + 
#   facet_grid(.~group, scales = "free") +
#   ylab("SSD time before index")+
#   xlab("Top 2 highest incidence Baddley et al.")+
#   ggtitle("SSD time before index by Baddley et al.")+ 
#   guides(fill="none")
# 
# 
# # SSD time west states
# 
# ssd_time <- read_csv(file = paste0(out_path,"ssd_duration_west_states.csv"))  
# 
# data <- ssd_time %>%
#   filter(!is.na(west_states)) %>% 
#   select(west_states, ssd_time_min, ssd_time_q1, ssd_time_median, ssd_time_q3, ssd_time_max) %>% 
#   pivot_longer(cols = ssd_time_min:ssd_time_max, values_to = "y", names_to = "group") %>% 
#   mutate(group = str_remove_all(group, "ssd_time_")) %>% 
#   mutate(new_lab=ifelse(west_states == 1, "yes", "no"))
# 
# data %>% 
#   mutate(group = factor(group, levels = c("min", "q1", "median", "q3", "max"))) %>% 
#   ggplot(aes(x = new_lab, y = y, fill = group))+
#   geom_bar(stat="identity", position=position_dodge()) + 
#   facet_grid(.~group, scales = "free") +
#   ylab("SSD time before index")+
#   xlab("West State")+
#   ggtitle("SSD time before index by West States vs. all others")+ 
#   guides(fill="none")
