
rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "tb"

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

# connect to database
# con <- DBI::dbConnect(RSQLite::SQLite(), paste0(delay_params$path,cond_name,".db"))

sim_in_path <- paste0(delay_params$path,"delay_results/")
sim_out_path <- paste0(delay_params$path,"delay_results/ssd_visit/")


if (!dir.exists(sim_out_path)) {
  dir.create(sim_out_path)
}

n_trials <- 1000

##############################
#### Simulation Functions ####
##############################

run_sim_draw_visits <- function(sim_data,overall_tm){
  tmp <- run_sim_draw_simple(sim_data)


  res1 <- distinct(tmp,obs)

  tmp2 <- tmp %>%
    group_by(patient_id) %>%
    summarise(duration=max(period),
              n_miss=n())

  count_duration <- function(val){
    tmp2 %>%
      summarise(n=sum(duration>=val))
  }

  res2 <- tmp2  %>%
    summarise(n_pat = n(),
              mean_dur = mean(duration),
              median_dur = median(duration),
              mean_n_miss = mean(n_miss),
              median_n_miss = median(n_miss))

  res3 <- tibble(duration_bin = sim_data$duration_bins) %>%
    mutate(res = map(duration_bin,count_duration)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp2),
           pct_all = 100*n/sim_data$total_patients)


  count_miss <- function(val){
    tmp2 %>%
      summarise(n = sum(n_miss>=val))
  }

  res4 <- tibble(miss_bin = sim_data$miss_bins) %>%
    mutate(res = map(miss_bin,count_miss)) %>%
    unnest(res) %>%
    mutate(pct_miss = 100*n/nrow(tmp2),
           pct_all = 100*n/sim_data$total_patients)

  return(list(obs = res1,
              stats = res2,
              dur_tab = res3,
              miss_tab = res4))
}

compute_export_stats <- function(sim_res,n_patients){
  ## Overall Stats ---------------------------------------------------------------

  sim_res_stats_means <- sim_res$sim_res_stats %>%
    summarise_all(funs(mean))

  sim_res_stats_upper <- sim_res$sim_res_stats %>%
    summarise_all(~quantile(.,probs = 0.975))

  sim_res_stats_lower <- sim_res$sim_res_stats %>%
    summarise_all(~quantile(.,probs = 0.025))

  sim_res_stats_means <- sim_res_stats_means %>%
    gather(key = key, value = mean) %>%
    inner_join(sim_res_stats_lower %>%
                 gather(key = key, value = low)) %>%
    inner_join(sim_res_stats_upper %>%
                 gather(key = key, value = high)) %>%
    filter(key != "trial") %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")"))

  # add in percent missed
  tmp <- sim_res_stats_means %>%
    filter(key == "n_pat") %>%
    mutate_at(vars(mean:high), ~100*./n_patients) %>%
    mutate(out = paste0(round(mean,2)," (",round(low,2),"-",round(high,2),")")) %>%
    mutate(key = "pct_miss")

  sim_res_stats_means <- bind_rows(tmp,sim_res_stats_means)

  rm(sim_res_stats_lower,sim_res_stats_upper,tmp)


  ## Duration stats --------------------------------------------------------------
  sim_res_dur_tab_means <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),funs(mean))


  sim_res_dur_tab_upper <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.975))

  sim_res_dur_tab_lower <- sim_res$sim_res_dur_tab %>%
    group_by(duration_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.025))


  sim_res_dur_tab_means <- sim_res_dur_tab_means %>%
    gather(key = key,value = mean, -duration_bin) %>%
    inner_join(sim_res_dur_tab_lower %>%
                 gather(key = key,value = low, -duration_bin)) %>%
    inner_join(sim_res_dur_tab_upper %>%
                 gather(key = key,value = high, -duration_bin)) %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")")) %>%
    select(duration_bin,key,out) %>%
    spread(key = key, value = out)

  rm(sim_res_dur_tab_lower,sim_res_dur_tab_upper)


  ## Miss stats ------------------------------------------------------------------
  sim_res_miss_tab_means <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),funs(mean))

  sim_res_miss_tab_upper <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.975))

  sim_res_miss_tab_lower <- sim_res$sim_res_miss_tab %>%
    group_by(miss_bin) %>%
    summarise_at(vars(n:pct_all),~quantile(.,probs = 0.025))

  sim_res_miss_tab_means <- sim_res_miss_tab_means %>%
    gather(key = key,value = mean, -miss_bin) %>%
    inner_join(sim_res_miss_tab_lower %>%
                 gather(key = key,value = low, -miss_bin)) %>%
    inner_join(sim_res_miss_tab_upper %>%
                 gather(key = key,value = high, -miss_bin)) %>%
    mutate(out = paste0(round(mean,3)," (",round(low,3),"-",round(high,3),")")) %>%
    select(miss_bin,key,out) %>%
    spread(key = key, value = out)

  rm(sim_res_miss_tab_lower,sim_res_miss_tab_upper)

  return(list(sim_res_stats = sim_res_stats_means,
              sim_res_dur_tab = sim_res_dur_tab_means,
              sim_res_miss_tab = sim_res_miss_tab_means))
}


generate_setting_counts <- function(obs_tm,sim_res){

  tmp <- obs_tm %>%
    inner_join(sim_res$sim_res_obs, by = "obs")

  tmp2 <- tmp %>%
    group_by(trial) %>%
    count(stdplac) %>%
    ungroup()

  tmp_fill <- tmp2 %>%
    distinct(stdplac) %>%
    mutate(trial = map(stdplac,~1:max(sim_res$sim_res_miss_tab$trial))) %>%
    unnest(trial)

  stdplac_res <- tmp2 %>%
    left_join(tmp_fill,.,by = c("stdplac", "trial")) %>%
    mutate(n=replace_na(n,0)) %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    group_by(stdplac) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975)) %>%
    mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
           pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    select(stdplac,n,pct,everything())


  tmp2 <- tmp %>%
    group_by(trial) %>%
    count(setting_type) %>%
    ungroup()

  setting_type_res <- tmp2 %>%
    group_by(trial) %>%
    mutate(pct = 100*n/sum(n)) %>%
    group_by(setting_type) %>%
    summarise(n_mean = mean(n),
              n_low = quantile(n,probs = 0.025),
              n_high = quantile(n,probs = 0.975),
              pct_mean = mean(pct),
              pct_low = quantile(pct,probs = 0.025),
              pct_high = quantile(pct,probs = 0.975)) %>%
    mutate(setting_type = setting_type_labels(setting_type)) %>%
    mutate(n = paste0(round(n_mean,2)," (",round(n_low,2),"-",round(n_high,2),")"),
           pct = paste0(round(pct_mean,2)," (",round(pct_low,2),"-",round(pct_high,2),")")) %>%
    select(setting_type,n,pct,everything())


  return(list(stdplac_res = stdplac_res,
              setting_type_res = setting_type_res))
}


run_sim <- function(sim_data,n_trials=1000, n){

  cluster <- parallel::makeCluster(35)

  parallel::clusterCall(cluster, function() library(tidyverse))
  parallel::clusterCall(cluster, function() library(delaySim))
  parallel::clusterExport(cluster,c("run_sim_draw_visits","sim_data"),
                          envir=environment())


  sim_res <- parallel::parLapply(cl = cluster,
                                 1:n_trials,
                                 function(x){run_sim_draw_visits(sim_data)})

  parallel::stopCluster(cluster)


  ## Aggregate sim results
  sim_res_obs <- map2(1:n_trials,sim_res,~mutate(.y$obs,trial=.x)) %>% bind_rows()
  sim_res_stats <- map2(1:n_trials,sim_res,~mutate(.y$stats,trial=.x)) %>% bind_rows()
  sim_res_dur_tab <- map2(1:n_trials,sim_res,~mutate(.y$dur_tab,trial=.x)) %>% bind_rows()
  sim_res_miss_tab <- map2(1:n_trials,sim_res,~mutate(.y$miss_tab,trial=.x)) %>% bind_rows()

  return(list(sim_res_obs = sim_res_obs,
              sim_res_stats = sim_res_stats,
              sim_res_dur_tab = sim_res_dur_tab,
              sim_res_miss_tab = sim_res_miss_tab))
}


######################
#### Prepare Data ####
######################

#### Load Data  ----------------------------------------------------------------

load(paste0(sim_in_path,"all_dx_visits.RData"))
load(paste0(sim_in_path,"delay_tm.RData"))

# load ssds
ssd_codes <- codeBuildr::load_ssd_codes(cond_name) %>%
  filter(type %in% c("icd9","icd10")) %>%
  mutate(dx=code,
         dx_ver = ifelse(type=="icd9",9L,10L)) %>%
  select(dx,dx_ver)

# index_dx_dates <- tbl(con,"index_dx_dates") %>% collect()

sim_tm <- all_dx_visits %>%
  mutate(period = -days_since_index) %>%
  inner_join(ssd_codes,by = c("dx", "dx_ver")) %>%
  distinct(enrolid,period,days_since_index) %>%
  inner_join(sim_obs,by = c("enrolid", "days_since_index")) %>% 
  rename(patient_id=enrolid) 

all_vis_count <- sim_tm %>%
  count(period) %>%
  filter(period>0)


# # merge observations into time map to extract visit types
obs_tm <- sim_tm %>%
  distinct(obs,days_since_index,enrolid=patient_id) %>%
  inner_join(tm,by = c("days_since_index", "enrolid")) %>%
  distinct(obs,stdplac,setting_type) %>%
  filter(setting_type!=4)

#### fit models ----------------------------------------------------------------
fit1 <- lm(n ~ period, filter(all_vis_count, period>delay_params$cp))
fit2 <- lm(n ~ log(period), filter(all_vis_count, period>delay_params$cp))
fit3 <- lm(n ~ poly(period,degree = 2), filter(all_vis_count, period>delay_params$cp))
fit4 <- lm(n ~ poly(period,degree = 3), filter(all_vis_count, period>delay_params$cp))

mse_res <- tibble(model_name = c("fit1","fit2","fit3","fit4"),
                  model = c("linear","exponential","quadratic","cubic"),
                  rmse = c(sqrt(sum(fit1$residuals^2)),
                           sqrt(sum(fit2$residuals^2)),
                           sqrt(sum(fit3$residuals^2)),
                           sqrt(sum(fit4$residuals^2)))) %>%
  mutate(label = paste("RMSE: ", round(rmse,2)))

y_pos <- .9*max(all_vis_count$n)
x_pos <- .8*delay_params$upper_bound

p1 <- all_vis_count %>%
  mutate(linear = predict(fit1,newdata = .)) %>%
  mutate(exponential = predict(fit2,newdata = .)) %>%
  mutate(quadratic = predict(fit3,newdata = .)) %>%
  mutate(cubic = predict(fit4,newdata = .)) %>%
  gather(key = model, value = value, -period, -n) %>%
  inner_join(mse_res) %>%
  ggplot(aes(period,n)) +
  geom_line() +
  scale_x_reverse() +
  geom_line(aes(y = value), color = "red") +
  facet_wrap(~model) +
  geom_vline(aes(xintercept = delay_params$cp), linetype =2) +
  theme_bw() +
  geom_text(data = mse_res,
            mapping = aes(x = x_pos, y = y_pos, label = label))


ggsave(filename = paste0(sim_out_path,"/expected_trends.pdf"),
       plot = p1)


#### Compute number of missed opportunities ------------------------------------

p2 <- all_vis_count %>%
  mutate(linear = predict(fit1,newdata = .)) %>%
  mutate(exponential = predict(fit2,newdata = .)) %>%
  mutate(quadratic = predict(fit3,newdata = .)) %>%
  mutate(cubic = predict(fit4,newdata = .)) %>%
  filter(period<=delay_params$cp) %>%
  gather(key = key,value = value, -period, -n) %>%
  mutate(num_miss = n-value) %>%
  mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Diagnosis") +
  theme_bw() +
  facet_wrap(~key) +
  ggtitle(paste0("Estimated Number of Missed Opportunities for ",cond_name))

ggsave(filename = paste0(sim_out_path,"/number_of_missed_opportunities.pdf"),
       plot = p2)


#################################
##### Evaluate Primary Model ####
#################################

# choose model
selected_model <- mse_res %>%
  filter(model!="exponential") %>%
  filter(rmse ==min(rmse))

fit <- eval(as.name(selected_model$model_name))


#### Compute number of missed opportunities ------------------------------------
sim_miss_bins <- all_vis_count %>%
  mutate(pred = predict(fit,newdata=.)) %>%
  mutate(num_miss = ifelse(n-pred>0,round(n-pred,0),0L)) %>%
  filter(period<=delay_params$cp) %>%
  select(period,num_miss)


#### Assemble final sim data ---------------------------------------------------

# assemble sim data
sim_data_all <- list(time_map = mutate(sim_tm,miss_ind=1),
                     miss_bins_visits = sim_miss_bins,
                     change_point = delay_params$cp,
                     total_patients = nrow(index_dx_dates),
                     miss_bins = delay_params$miss_bins,
                     duration_bins = delay_params$duration_bins)


#### Run Simulation ------------------------------------------------------------
sim_res_main <-  run_sim(sim_data_all,n_trials = n_trials)

export_stats <- compute_export_stats(sim_res = sim_res_main,
                                     n_patients = sim_data_all$total_patients)

setting_counts <- generate_setting_counts(obs_tm = obs_tm,
                                          sim_res = sim_res_main)

# Save output
optimal_model_res <- c(export_stats,setting_counts)

save(optimal_model_res,file = paste0(sim_out_path,"/optimal_model_res.RData"))

# pull out distinct observations drawn
sim_res <- sim_res_main$sim_res_obs

# Save distinct observations in simulation
sim_res_sim_obs <- sim_res %>% 
  distinct(obs) %>% 
  inner_join(sim_obs) 

# Save simulation draws for later analysis
save(sim_res,sim_res_sim_obs,file = paste0(sim_out_path,"/sim_res.RData"))


################################
##### Evaluate Other Models ####
################################

other_mods <- mse_res %>%
  anti_join(selected_model)

out <- vector("list",3)

names(out) <- other_mods$model

for (i in 1:nrow(other_mods)){

  print(other_mods$model[i])

  res <- list()

  # extract model
  fit <- eval(as.name(other_mods$model_name[i]))

  # compute miss bins
  sim_miss_bins_other <- all_vis_count %>%
    mutate(pred = predict(fit,newdata=.)) %>%
    mutate(num_miss = ifelse(n-pred>0,round(n-pred,0),0L)) %>%
    filter(period<=delay_params$cp) %>%
    select(period,num_miss)

  res[["sim_miss_bins"]] <- sim_miss_bins_other
  res[["total_miss"]] <- sum(sim_miss_bins_other$num_miss)


  sim_data_other <- list(time_map = mutate(sim_tm,miss_ind=1),
                         miss_bins_visits = sim_miss_bins_other,
                         change_point = delay_params$cp,
                         total_patients = nrow(index_dx_dates),
                         miss_bins = delay_params$miss_bins,
                         duration_bins = delay_params$duration_bins)


  #### Run Simulation ------------------------------------------------------------
  sim_res_other <-  run_sim(sim_data_other)

  export_stats <- compute_export_stats(sim_res = sim_res_other,
                                       n_patients = sim_data_other$total_patients)

  setting_counts <- generate_setting_counts(obs_tm = obs_tm,
                                            sim_res = sim_res_other)

  res[["model_res"]] <-  c(export_stats,setting_counts)

  out[[i]] <- res


}

# Save Results
other_model_res <- out

save(other_model_res,file = paste0(sim_out_path,"/other_model_res.RData"))


######################################
#### Evaluate Other Change-points ####
######################################



##########################################
#### Extract simulation probabilities ####
##########################################

# sim_res <- sim_res_main$sim_res_obs
# 
# # visit probability
# vis_miss_probs <- sim_res %>%
#   count(obs) %>%
#   arrange(desc(n)) %>%
#   mutate(prob = n/n_trials)
# 
# vis_miss_probs <- vis_miss_probs %>%
#   inner_join(sim_tm %>%
#                distinct(obs,period,patient_id))
# 
# # patient probability
# pat_miss_probs <- sim_res %>%
#   inner_join(sim_tm %>%
#                distinct(obs,period,patient_id)) %>%
#   distinct(trial,patient_id) %>%
#   count(patient_id) %>%
#   arrange(desc(n)) %>%
#   mutate(miss_prob = n/n_trials)
# 
# # patient duration probability
# pat_dur_probs <- sim_res %>%
#   inner_join(sim_tm %>%
#                distinct(obs,period,patient_id)) %>%
#   group_by(patient_id,trial) %>%
#   summarise(duration = max(period)) %>%
#   count(patient_id,duration) %>%
#   ungroup() %>%
#   mutate(prob = n/n_trials)
# 
# sim_res_probs <- list(vis_miss_probs = vis_miss_probs,
#                       pat_miss_probs = pat_miss_probs,
#                       pat_dur_probs = pat_dur_probs)
# 
# save(sim_res_probs, file = paste0(sim_out_path,"/sim_res_probs.RData"))


############################################
#### Create Readme of analysis overview ####
############################################

# add readme noting build info
heading <- paste0("Simulation for ",cond_name," using any visits")
sim_date <- paste0("Simulation Date: ",Sys.time())
trial_size <- paste0("Simulation Trials: ",n_trials)
population_size <- paste0("Total Population: ",n_trials)
change_point_used <- paste0("Change-point used: ",delay_params$cp," days")
optimal_model <- paste0("Selected Model: ",selected_model$model)
tmp <- sim_data_all$miss_bins_visits %>% summarise(n_miss =sum(num_miss)) %>% .$n_miss
total_miss <- paste0("Total Number of misses: ",tmp)
tmp1 <- names(other_model_res)
tmp2 <- map(other_model_res,~.$total_miss)
total_miss_other <- c("Total miss other models: \n",paste0("  ",tmp1," - ",tmp2," \n"))

# collect results
tmp <- map(other_model_res,~.$model_res$sim_res_stats)
tmp[[selected_model$model]] <- optimal_model_res$sim_res_stats

pct_miss_summary <- c("Percent Miss: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[1])," \n"))
mean_dur_summary <- c("Mean Duration: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[3])," \n"))
median_dur_summary <- c("Median Duration: \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[4])," \n"))
mean_n_miss_summary <- c("Mean Number of Misses (per patient): \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[5])," \n"))
median_n_miss_summary <- c("Median Number of Misses (per patient): \n",paste0("  ",names(tmp)," - ",map(tmp,~.$out[6])," \n"))


fileConn<-file(paste0(sim_out_path,"readme.txt"))
writeLines(c(heading,
             sim_date,
             trial_size,
             population_size,
             change_point_used,
             optimal_model,
             total_miss,
             total_miss_other,
             pct_miss_summary,
             mean_dur_summary,
             median_dur_summary,
             mean_n_miss_summary,
             median_n_miss_summary), fileConn)
close(fileConn)





