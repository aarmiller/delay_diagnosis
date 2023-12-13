library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)
# library(changepoint)
# library(trend)

params <- list()
params$cond <- "dengue"

# load delay_parms
load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[params$cond]]

cond_name <- filter(codeBuildr::avail_disease_codes(),name==params$cond)$description

in_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",params$cond,"/delay_results")
# in_path <- paste0("/Volumes/Statepi_Diagnosis/prelim_results/",params$cond,"/delay_results")

load(paste0(in_path,"/all_dx_visits.RData"))

# load SSD codes
ssd_codes <- bind_rows(codeBuildr::load_ssd_codes(params$cond) %>%
                         filter(type == "icd9") %>%
                         select(dx=code) %>%
                         mutate(dx_ver=9L),
                       codeBuildr::load_ssd_codes(params$cond) %>%
                         filter(type == "icd10") %>%
                         select(dx=code) %>%
                         mutate(dx_ver=10L))


# range of CP to evaluate
cp_range <- 28:(round(delay_params$upper_bound*(2/3)))


#### Functions -----------------------------------------------------------------

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

fit_range <- function(data,cp_range){
  tibble(lb = cp_range) %>% 
    mutate(res = map(lb,~fit_trends(data,.)$cp_res)) %>% 
    unnest(res)
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

###################
#### Prep Data ####
###################

# Note this step is redundant for now...this is where the final cohort would be filtered
# patient_ids <- index_dx_dates %>% distinct(patient_id)
# all_dx_visits <- all_dx_visits %>% inner_join(patient_ids)
# 
# visit_counts <- all_dx_visits %>%
#   distinct(patient_id,dx_ver,days_since_index) %>%
#   count(dx_ver,days_since_index)


count_data_all <- visit_counts %>% 
  filter(is.na(dx_ver)) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  mutate(dow = paste0("dow_",dow))

# populate missing valules in visit counts (i.e., assign 0 to days missing)
count_data_ssd <- all_dx_visits %>%
  inner_join(ssd_codes, by = join_by(dx, dx_ver)) %>%
  distinct(patient_id,days_since_index) %>%
  count(days_since_index) %>%
  left_join(tibble(days_since_index=min(visit_counts$days_since_index):-1),., by = "days_since_index") %>% # in case there are 0 days
  mutate(n = replace_na(n,0L)) %>% 
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index) %>% 
  filter(days_since_index<0) %>% 
  mutate(period = -days_since_index) %>% 
  mutate(dow = as.factor(period %% 7)) %>% 
  mutate(dow = paste0("dow_",dow))


### Prep Bootstrap Data --------------------------------------------------------

# sample overall set of patients for each bootstrap
tmp <- tibble(boot_trial=1:100) %>% 
  mutate(boot_sample = map(boot_trial,~tibble(patient_id = sample(index_dx_dates$patient_id, replace = TRUE)))) %>% 
  mutate(boot_sample = map(boot_sample, ~mutate(.,boot_id = row_number())))

# tmp$boot_sample[[1]]

#### compute counts for all visits ---------------------------------------------
# For each set of patients selected pull the time map visit counts
tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "patient_id", relationship = "many-to-many") %>%
                        mutate(period = -days_since_index) %>%
                        distinct(patient_id,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("patient_id", "days_since_index"))))

tmp <- tmp %>% 
  mutate(all_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7)) %>% 
                               mutate(dow = paste0("dow_",dow))))

#### compute counts for ssd visits ---------------------------------------------
# for each set of patients pull the ssd counts

tmp <- tmp %>% 
  mutate(sim_tm = map(boot_sample,
                      ~inner_join(.,all_dx_visits, by = "patient_id", relationship = "many-to-many") %>%
                        mutate(period = -days_since_index) %>%
                        inner_join(ssd_codes,by = c("dx", "dx_ver")) %>%
                        distinct(patient_id,period,days_since_index,boot_id) %>%
                        inner_join(sim_obs,by = c("patient_id", "days_since_index"))))

tmp <- tmp %>% 
  mutate(ssd_vis_count = map(sim_tm,
                             ~count(.,period) %>%
                               filter(period>0) %>% 
                               mutate(dow = as.factor(period %% 7)) %>% 
                               mutate(dow = paste0("dow_",dow))))

# remove the sim_tm data (no longer needed)
tmp <- tmp %>% 
  select(boot_trial,boot_sample,all_vis_count,ssd_vis_count)

boot_counts <- tmp
rm(tmp)



###################
### Analyze CP ####
###################


boot_counts %>% 
  select(boot_trial,ssd_vis_count) %>% 
  unnest(ssd_vis_count) %>% 
  inner_join(select(count_data_ssd,period,n_ssd=n)) %>% 
  ggplot(aes(period,n,group = boot_trial)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(y = n_ssd), size = 0.3, color = "red") +
  theme_minimal() +
  scale_x_reverse()

### Analyze Change-points ------------------------------------------------------


tmp <- fit_range(boot_counts$ssd_vis_count[[1]],28:(round(delay_params$upper_bound*(2/3))))


ssd_boot_res <- boot_counts %>% 
  select(boot_trial,ssd_vis_count) %>% 
  mutate(cp_res = map(ssd_vis_count,~fit_range(.,cp_range)))


# Save output for transfer
save(boot_counts,count_data_all,count_data_ssd,ssd_boot_res,file = "/Shared/Statepi_Diagnosis/prelim_results/dengue/change_point_results/new_cp_res.RData")

load("/Volumes/Statepi_Diagnosis/prelim_results/dengue/change_point_results/new_cp_res.RData")

ssd_boot_res$cp_res

ssd_boot_res %>% 
  select(boot_trial,ssd_vis_count) %>% 
  unnest(ssd_vis_count) %>% 
  ggplot(aes(period,n,group = boot_trial)) +
  geom_line(alpha = 0.2) +
  theme_minimal() +
  scale_x_reverse()

# Most common implied CP
tmp <- ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  filter(n == max(n)) %>% 
  ungroup()
tmp


# Distribution of implied CPs
ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>% 
  inner_join(select(tmp,model,mc_implied_cp = implied_cp)) %>% 
  ggplot(aes(implied_cp)) +
  geom_histogram() +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = mc_implied_cp), color = "red")

# average deviation from lower bound for the most common implied cp
select(tmp,model,implied_cp) %>% 
  inner_join(ssd_boot_res %>% 
               select(boot_trial,cp_res) %>% 
               unnest(cp_res)) %>% 
  mutate(dev = lb-implied_cp) %>% 
  group_by(model) %>% 
  summarise(mean_dev = mean(dev), median_dev = median(dev),
            mean_abs_dev = mean(abs(dev))) %>% 
  left_join(tmp,.) %>% 
  select(-n)

# Now restrict to change-points within some abs deviation

# 5 days
ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>%
  filter(abs(lb-implied_cp)<=5) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  filter(n == max(n)) %>% 
  ungroup()

# 1 day
ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>%
  filter(abs(lb-implied_cp)<=1) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  filter(n == max(n)) %>% 
  ungroup()

ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>%
  filter(abs(lb-implied_cp)<=5) %>% 
  ggplot(aes(implied_cp)) +
  geom_histogram() +
  facet_wrap(~model)


ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>%
  filter(abs(lb-implied_cp)<=1) %>% 
  ggplot(aes(implied_cp)) +
  geom_histogram() +
  facet_wrap(~model)

# Reapply to final model
fit_trends(count_data_ssd,30)
fit_trends(count_data_ssd,18)

fit_trends(count_data_ssd,30) %>% 
  .$pred %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y=pred)) +
  facet_wrap(~model) 
  


ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>%
  filter(abs(lb-implied_cp)<=5) %>% 
  filter(implied_cp == 30) %>% 
  group_by(model) %>% 
  summarise(mean_mse = mean(mse))

tmp <- ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n))

ssd_boot_res$cp_res[[1]] %>% 
  mutate(cp_dev = abs(lb-implied_cp)) %>% 
  filter(cp_dev<=5) %>% 
  count(implied_cp)





ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>% 
  mutate(cp_dev = abs(lb-implied_cp)) %>% 
  filter(cp_dev<=5) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  ggplot(aes(implied_cp,n)) +
  geom_histogram(stat = "identity") +
  facet_wrap(~model, scales = "free")

ssd_boot_res %>% 
  select(boot_trial,cp_res) %>% 
  unnest(cp_res) %>% 
  mutate(cp_dev = abs(lb-implied_cp)) %>% 
  filter(cp_dev<=0) %>% 
  count(model,implied_cp) %>% 
  group_by(model) %>% 
  mutate(frac = 100*n/sum(n)) %>% 
  ggplot(aes(implied_cp,n)) +
  geom_histogram(stat = "identity") +
  facet_wrap(~model, scales = "free")



tibble(lb = 28:(round(delay_params$upper_bound*(2/3)))) %>% 
  mutate(res = map(lb,~fit_trends(count_data_ssd,.)$cp_res))
  


fit_trends(count_data_ssd,128)$cp_res


tmp <- fit_trends(count_data_ssd,128)

round(delay_params$upper_bound*(2/3))

fits_ssd <- tibble(cp = 28:(round(delay_params$upper_bound*(2/3)))) %>% 
  mutate(res = map(cp,~fit_trends(count_data_ssd,.)))

fits_ssd <- fits_ssd %>% 
  mutate(cp_res = map(res,~.$cp_res),
         pred = map(res,~.$pred)) %>% 
  select(-res)


fit_select <- tibble(cp_deviation = 0:15) %>% 
  mutate(cp_fp=map(cp_deviation,~0:20)) %>% 
  unnest(cp_fp) %>% 
  mutate(res = map2(cp_deviation,cp_fp,~find_cp(fits_ssd,.x,.y)))

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


tmp$pred %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  scale_x_reverse() +
  theme_minimal() +
  geom_vline(aes(xintercept = implied_cp), linetype = 2) +
  xlab("Days before Diagnosis") +
  ylab("Number of Visits") +
  facet_wrap(~model)

###############
#### Plots ####
###############

### Comparison of All and SSD visit trends

bind_rows(mutate(count_data_all, group = "all"),
          mutate(count_data_ssd, group = "ssd")) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_wrap(~group, scales = "free_y") 
