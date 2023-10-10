
rm(list = ls())
sarcoid_lung

#### Load Data

sarcoid_lung <- read_csv("cp_approaches/data/sarcoid_lung_ssd_counts.csv") 
sarcoid_skin <- read_csv("cp_approaches/data/sarcoid_skin_ssd_counts.csv")

load("/Volumes/AML/tmp_transfer/saroid_lung_bootstrap.RData")
sarcoid_lung_bootstrap <- tmp
rm(tmp)

sarcoid_lung_bootstrap %>% 
  select(boot_trial,ssd_vis_count) %>% 
  unnest(ssd_vis_count) %>% 
  mutate(boot_trial=as.factor(boot_trial)) %>% 
  inner_join(select(sarcoid_lung,period,n_aggregate=n)) %>% 
  ggplot(aes(period,n,group = boot_trial)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  scale_x_reverse() 
  # geom_line(aes(y = n_aggregate), color = "red", size = .1)


#### Functions -------------------------------------

fit_trends <- function(count_data,lower_bound){
  
  # filter to training data
  tmp_data <- count_data %>% 
    filter(period>=lower_bound)
  
  # fit models
  fit_lm <- lm(n~period+dow, data = tmp_data)
  fit_quad <- lm(n~poly(period,2)+dow, data = tmp_data)
  fit_cube <- lm(n~poly(period,3)+dow, data = tmp_data)
  fit_exp <- lm(n~log(period)+dow, data = tmp_data)
  
  # get fitted values
  pred_data_tmp <- count_data %>% 
    mutate(lm = predict(fit_lm,newdata =.),
           quad = predict(fit_quad,newdata =.),
           cube = predict(fit_cube,newdata =.),
           exp = predict(fit_exp,newdata =.)) 
  
  find_pred_bound_cp <- function(pred_data){
    tmp_lm <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>lm)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_quad <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>quad)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_cube <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>cube)==row_number()) %>% 
      filter(above_cp) 
    
    tmp_exp <- pred_data %>% 
      arrange(period) %>% 
      mutate(above_cp = cumsum(n>exp)==row_number()) %>% 
      filter(above_cp) 
    
    tibble(model = c("lm","quad","cube","exp"),
           implied_cp = c(ifelse(nrow(tmp_lm)==0,NA,max(tmp_lm$period)),
                          ifelse(nrow(tmp_quad)==0,NA,max(tmp_quad$period)),
                          ifelse(nrow(tmp_cube)==0,NA,max(tmp_cube$period)),
                          ifelse(nrow(tmp_exp)==0,NA,max(tmp_exp$period))))
  }
  
  implied_cps <- find_pred_bound_cp(pred_data_tmp)
  
  pred_data_out <- pred_data_tmp %>% 
    select(period,n,lm:exp) %>% 
    gather(key = "model", value = "pred", -period, -n) %>% 
    inner_join(implied_cps, by = "model") 
  
  
  find_fit_measures <- function(pred_data,bound_val){
    pred_data %>%
      group_by(model) %>%
      filter(period>=implied_cp) %>% 
      summarise(mse=mean((n-pred)^2))
  }
  
  out1 <- inner_join(implied_cps,find_fit_measures(pred_data_out,lower_bound), by = "model")
  
  return(list(cp_res=out1,
              pred = pred_data_out))
  
}

find_cp <- function(data,cp_deviation,cp_fp){
  
  tmp <- data %>% 
    select(cp,cp_res) %>% 
    unnest(cp_res) %>% 
    mutate(cp_diff = abs(cp-implied_cp)) %>% 
    filter(implied_cp<=(cp-cp_deviation),
           cp_diff<=cp_fp)            # allow implied change-point to deviate to the left of the provided cp
  
  if (nrow(tmp)>0){ 
    return(group_by(tmp,model) %>% 
             filter(mse==min(mse)))
  } else {
    tmp
  }
}


#### Evaluate the raw data ---------------------------

tmp <- fit_trends(sarcoid_lung,128)

bind_rows(sarcoid_lung %>% 
            mutate(series = "Lung"),
          sarcoid_skin %>% 
            mutate(series = "Skin")) %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  facet_wrap(~series, scales = "free_y") +
  theme_minimal() +
  scale_x_reverse() +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") +
  theme(strip.text = element_text(size = 12))

tmp$pred %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  # geom_line(aes(y = pred, color = model)) +
  scale_x_reverse() +
  theme_minimal() +
  geom_vline(aes(xintercept = 128), linetype = 2) +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") 

tmp <- tibble(cp = 100:500) %>% 
  mutate(res = map(cp,~fit_trends(sarcoid_lung,.)))

tmp <- tmp %>% 
  mutate(cp_res = map(res,~.$cp_res),
         pred = map(res,~.$pred)) %>% 
  select(-res)


tmp %>% 
  select(cp,cp_res) %>% 
  unnest(cp_res) %>% 
  filter(model=="exp") %>% 
  filter(cp==157)

fit_trends(sarcoid_lung,300)$pred %>% 
  filter(model == "exp") %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred), color = "green") +
  scale_x_reverse() +
  theme_minimal() +
  geom_vline(aes(xintercept = 300), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 183), linetype = 2, color = "red") +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") 

tmp %>% 
  select(cp,cp_res) %>% 
  unnest(cp_res) %>% 
  filter(model=="exp",
         cp == 157)
  
fit_trends(sarcoid_lung,183)$pred %>% 
  filter(model == "exp")

fit_trends(sarcoid_lung,183)$pred %>% 
  filter(model == "exp") %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred), color = "green") +
  scale_x_reverse() +
  theme_minimal() +
  geom_vline(aes(xintercept = 183), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 157), linetype = 2, color = "red") +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") 

fit_trends(sarcoid_lung,157)$pred %>% 
  filter(model == "exp")

fit_trends(sarcoid_lung,157)$pred %>% 
  filter(model == "exp") %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  # geom_line(aes(y = pred), color = "green") +
  scale_x_reverse() +
  theme_minimal() +
  # geom_vline(aes(xintercept = 157), linetype = 2, color = "blue") +
  # geom_vline(aes(xintercept = 135), linetype = 2, color = "red") +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") 



fit_trends(sarcoid_lung,135)$pred %>% 
  filter(model == "exp")

tmp %>% 
  select(cp,cp_res) %>% 
  unnest(cp_res) %>% 
  filter(model=="exp",
         cp == 135)

fit_trends(sarcoid_lung,135)$pred %>% 
  filter(model == "exp") %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred), color = "green") +
  scale_x_reverse() +
  theme_minimal() +
  geom_vline(aes(xintercept = 135), linetype = 2, color = "blue") +
  geom_vline(aes(xintercept = 128), linetype = 2, color = "red") +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") 



tmp_res <- tibble(cp_deviation = 0:15) %>% 
  mutate(cp_fp=map(cp_deviation,~0:20)) %>% 
  unnest(cp_fp) %>% 
  mutate(res = map2(cp_deviation,cp_fp,~find_cp(tmp,.x,.y)))

tmp_res %>% 
  unnest(res)  %>% 
  inner_join(tibble(model = c("cube","exp","lm","quad"),
                    model_label = c("Cubic","Exponential","Linnear","Quadratic"))) %>% 
  ggplot(aes(cp_deviation,cp_fp,fill = mse)) +
  geom_tile() +
  geom_text(aes(label = implied_cp)) +
  scale_fill_gradient(low = "red", high = "white") +
  facet_wrap(~model_label) +
  theme_minimal() +
  ylab("Fixed Point Deviation") +
  xlab("Left Side Deviation") +
  theme(strip.text = element_text(size = 12))


fin_mods <- tmp_res %>% 
  unnest(res) %>% 
  group_by(model) %>% 
  filter(cp_deviation<10) %>% 
  filter(cp_fp<10) %>% 
  filter(mse == min(mse)) %>% 
  filter(cp_deviation==min(cp_deviation)) %>% 
  filter(cp_fp==min(cp_fp))

fin_mods <- fin_mods %>%
  inner_join(tibble(model = c("cube","exp","lm","quad"),
                    model_label = c("Cubic","Exponential","Linnear","Quadratic"))) %>% 
  mutate(model_label = paste0(model_label,", CP=",implied_cp,", MSE=",round(mse,2)))
  

tmp %>% 
  select(cp,pred) %>% 
  unnest(pred) %>% 
  inner_join(fin_mods) %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred), color = "red") +
  facet_wrap(~model_label) +
  geom_vline(aes(xintercept=implied_cp), linetype = 2) +
  theme_minimal() +
  scale_x_reverse() +
  theme(strip.text = element_text(size = 12)) +
  ylab("Number of Visits") +
  xlab("Days Before Diagnosis")

fin_mods

11111tmp_res %>% 
  unnest(res) %>% 
  mutate()

tmp_res %>% 
  unnest(res) %>% 
  filter(model=="pred_exp")



tmp %>% 
  select(cp,cp_res) %>% 
  unnest(cp_res) %>% 
  mutate(cp_diff = abs(cp-implied_cp)) %>% 
  filter(implied_cp<=(cp-1)) %>%            # allow implied change-point to deviate to the left of the provided cp
  filter(cp_diff<=1) %>%                   # allow the cp and implied cp to differ by given amount
  group_by(model) %>% 
  filter(mse==min(mse))
arrange(mse)



tmp2 <- fit_trends(count_data_ssd,383) 

tmp2$pred %>% 
  select(period,n,pred_lm:pred_exp) %>% 
  gather(key = "model", value = "pred", -period, -n) %>%
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = pred, color = model))



#### Evaluate stability over raw counts ----------------------------------------

fit_periods <- tibble(cp = 100:500) %>% 
  mutate(res = map(cp,~fit_trends(sarcoid_lung,.))) %>% 
  mutate(cp_res = map(res,~.$cp_res),
         pred = map(res,~.$pred)) %>% 
  select(-res)

tmp_res <- tibble(cp_deviation = 0:15) %>% 
  mutate(cp_fp=map(cp_deviation,~0:20)) %>% 
  unnest(cp_fp) %>% 
  mutate(res = map2(cp_deviation,cp_fp,~find_cp(fit_periods,.x,.y)))

tmp_res %>% 
  unnest(res)  %>% 
  inner_join(tibble(model = c("cube","exp","lm","quad"),
                    model_label = c("Cubic","Exponential","Linnear","Quadratic"))) %>% 
  ggplot(aes(cp_deviation,cp_fp,fill = mse)) +
  geom_tile() +
  geom_text(aes(label = implied_cp)) +
  scale_fill_gradient(low = "red", high = "white") +
  facet_wrap(~model_label) +
  theme_minimal() +
  ylab("Fixed Point Deviation") +
  xlab("Left Side Deviation") +
  theme(strip.text = element_text(size = 12))

tmp_res %>% 
  unnest(res) %>% 
  filter(cp_deviation==0) %>% 
  group_by(model) %>% 
  filter(mse==min(mse)) %>% 
  ungroup() %>% 
  distinct()


fin_mods <- tmp_res %>% 
  unnest(res) %>% 
  group_by(model) %>% 
  filter(cp == implied_cp) %>% 
  filter(cp_deviation<10) %>% 
  filter(cp_fp<10) %>% 
  distinct(cp,model,implied_cp) %>% 
  ungroup()

fin_mods %>% count(implied_cp,model)


fin_mods %>% inner_join(tmp_res %>% 
                          unnest(res)) %>% 
  group_by(model) %>% 
  filter(mse==min(mse)) %>% 
  ungroup() %>% 
  distinct(cp,model,mse)

#### Evaluate stability over bootstraps ----------------------------------------

cps_to_eval <- distinct(fin_mods,cp) %>% .$cp

tmp <- fit_trends(sarcoid_lung,128)

tmp <- sarcoid_lung_bootstrap %>% 
  select(boot_trial,ssd_vis_count) %>% 
  mutate(cp = map(boot_trial, ~cps_to_eval)) %>% 
  unnest(cp) %>% 
  mutate(res = map2(ssd_vis_count,cp,~fit_trends(.x,.y)))

tmp2 <- tmp %>% 
  mutate(cp_res = map(res,~.$cp_res)) %>% 
  select(boot_trial,cp,cp_res) %>% 
  unnest(cp_res)

tmp3 <- rename(tmp2,implied_cp_boot=implied_cp) %>% 
  inner_join(fin_mods)

tmp4 <- tmp3 %>%
  group_by(model,implied_cp_boot) %>% 
  count(model,implied_cp_boot) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  arrange(model,desc(frac))

tmp4 %>% 
  filter(model == "exp")

tmp3 %>% 
  mutate(cp_dev = cp-implied_cp_boot) %>% 
  ggplot(aes(cp_dev)) +
  geom_histogram() +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal()


tmp3 %>% 
  mutate(cp_dev = cp-implied_cp_boot) %>% 
  filter(cp == implied_cp) %>% 
  ggplot(aes(cp_dev)) +
  geom_histogram() +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal()

tmp3 %>% 
  filter(implied_cp == cp) %>% 
  mutate(cp_dev = abs(cp-implied_cp_boot)) %>% 
  group_by(model) %>% 
  summarise(`Mean Dev` = mean(cp_dev),
            `Median Dev.` = median(cp_dev),
            `Mean MSE` = mean(mse))


eval_boot_mods <- function(data){
  tmp1 <- data %>%
    group_by(model,implied_cp_boot) %>% 
    summarise(n_trials = n(),
              mean_mse = mean(mse)) %>% 
    ungroup()
  
  tmp2 <- data %>% 
    filter(cp==implied_cp) %>% 
    group_by(model) %>% 
    summarise(mse_fp = mean(mse)) %>% 
    ungroup()
  
  return(list(trial_results = tmp1,
              fp_results))
}

eval_boot_mods(tmp3)

tmp3 %>% 
  group_by(model,implied_cp) %>% 
  summarise(mean_mse = mean(mse))

tmp3 %>% 
  filter(model == "cube") %>% 
  arrange(mse) %>% 
  summarise(mean)

tmp3 %>% 
  filter(cp==implied_cp) %>% 
  group_by(model) %>% 
  summarise(mse = mean(mse))

fit_trends(sarcoid_lung,128)
