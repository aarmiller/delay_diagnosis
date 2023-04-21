
############################################################
#### Script to Find Change-Points and Related Quantities ###
############################################################

### Load packages and define necessary quantities ###

rm(list=ls())
library(tidyverse)
# library(icdplus)
# library(RSQLite)
library(bit64)
library(lubridate)
# library(stringr)
# library(parallel)
# install.packages("reReg")
# library(reReg)
library(delaySim)
library(smallDB)
# library(DBI)
library(forecast)
library(predictBoundCP)
library(knitr)
library(rmarkdown)
library(kableExtra)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
codeBuildr::avail_ssd_codes()
library(codeBuildr)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "stroke"

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

in_path <- paste0(delay_params$path,"delay_results")
out_path <- paste0(delay_params$path,"change_point_results/")

if (!dir.exists(out_path)){
  dir.create(out_path)
}

##  Set days before to count delay ##
days_before <- 1L

##  Set by_days to count delay ##
by_days <- 1L

########################
### Load in all data ###
########################

load(paste0(in_path,"/all_dx_visits.RData"))

#Change this to read from wherever the data is stored for you
visit_counts <- visit_counts %>%
  filter(is.na(dx_ver)) %>%
  filter(days_since_index<0) %>%
  mutate(period = -days_since_index) %>%
  select(period,n,days_since_index)

# visit_counts %>% arrange(period)

#Get max days we are looking before in this dataset
#Can also set this manually if you want

max_days <- -1*min(visit_counts$days_since_index)

##################################################
### Calculate data for all changepoint methods ###
##################################################

all_methods <- tibble(method = c("lm", "lm_quad", "lm_cube", "quad", "cube", "exp",
                                 # "spline",
                                 "pettitt", "cusum"
                                 # "MSE", "RMSE", "MAE", "MSLE", "RMSLE"
                                 )) %>%
  mutate(week_period = map(method, ~c(TRUE, FALSE))) %>% unnest(cols = week_period) %>%
  mutate(auto_reg = map(method, ~c(FALSE))) %>% unnest(cols = auto_reg) %>%
  # mutate(auto_reg = map(method, ~c(TRUE, FALSE))) %>% unnest(cols = auto_reg) %>%
  mutate(var = map(method, ~c("n"))) %>% unnest(cols = var)


final_prior_visit_counts <- all_methods %>% mutate(counts = pmap(list(method),~filter(visit_counts,period<=max_days)))



###Fitting models###
#run set change point function
results <- final_prior_visit_counts %>%
  mutate(out=pmap(list(counts, method, week_period, auto_reg),
                  ~delaySim::find_change_point(data = ..1, method = ..2,
                                               var_name = "n",
                                               week_period = ..3, specify_cp = NULL,
                                               auto_reg = ..4)))

results <- results %>% mutate(expected_trend=map(out,~.$pred %>% select(period, t, pred1)))

#extract change points
final_results <- results %>%
  mutate(change_point=map(out, ~.$change_point$period)) %>%
  unnest(change_point)

###Aaron, I leave this in here if you want these plots, since they are
###already calculated by the find_change_point function automatically
#Uncomment to include plots in the final_results object
#create plots of change points
#final_results <- final_results %>%
#  mutate(graphs=map(out, ~.$cp_plot))

final_results <- final_results %>% select(-counts, -out)


#################We now do a similar procedure, but for the prediction bound method#############################




#Create a function to generate this data
generate_pred_bound_results <- function(data, week_period, auto_reg , set_low){

  if(week_period | auto_reg){
    this_model <- "lm_period"
  } else{
    this_model <- "lm"
  }

  ## Prediction bound method with 90% bounds ##
  pred_boud_lm <- predictBoundCP::fit_cp_model(count_data = data %>%
                                                 mutate(days_since_dx = period * -1) %>%
                                                 select(days_since_dx, n),
                                               low = set_low, high = max_days,
                                               model = this_model,
                                               level = 0.90,
                                               plot = FALSE)

  #run set change point function to generate the results we seek
  results <- data %>%
    delaySim::find_change_point(method = "cusum",
                                specify_cp = pred_boud_lm$cp,
                                var_name = "n",
                                week_period = week_period,
                                auto_reg = auto_reg)

  #Change plots labels
  results$cp_plot$labels$title <- paste0("Method = pred_bound, week_effect = ",week_period,", auto_reg=",auto_reg)

  return(results)
}

###Now, we prepare data and run this function
# specify low to be 10% greater than cp found by cusum
#run set pred bound function
results_pred_bound <- final_prior_visit_counts %>% filter(method == "cusum") %>%
  mutate(low = pmap(list(week_period, auto_reg, var), ~round(final_results %>%
                                                               filter(method == "cusum",
                                                                      week_period == ..1,
                                                                      auto_reg == ..2,
                                                                      var == ..3) %>% .$change_point * 1.10, 0)),
         out = pmap(list(counts, week_period, auto_reg, low),
                    ~generate_pred_bound_results(data = ..1,
                                                 week_period = ..2,
                                                 auto_reg = ..3,
                                                 set_low = ..4)))

results <- results_pred_bound %>% mutate(expected_trend=map(out,~.$pred %>% select(period, t, pred1)))

#extract change points
final_results2 <- results %>%
  mutate(change_point=map(out, ~.$change_point$period)) %>%
  unnest(change_point)


###Aaron, I leave this in here if you want these plots, since they are
###already calculated by the find_change_point function automatically
#Uncomment to include plots in the final_results object
#create plots of change points
#final_results2 <- final_results2 %>%
#  mutate(graphs=map(out, ~.$cp_plot))

final_results2 <- final_results2 %>% select(-out, -counts, -low) %>% mutate(method = "pred_bound")

########We append this data to our other datasets##########

final_results <- rbind(final_results, final_results2)

#Save these results
second_pass_pred_bound <- function(pred_data,count_data){
  tmp <- pred_data %>%
    inner_join(count_data,by = "period") %>%
    arrange(period) %>%
    mutate(above_cp=cumsum(n>pred1)==row_number()) %>%
    filter(above_cp) %>%
    summarise(cp=max(period))

  as.integer(tmp$cp)
}


final_results <- final_results %>%
  mutate(cp_pred_bound = map_int(expected_trend,~second_pass_pred_bound(.,visit_counts)))

#### Export Final change-points ------------------------------------------------
final_results %>%
  select(method,periodicity=week_period,auto_reg,change_point,cp_pred_bound) %>%
  write_csv(paste0(out_path,"any_visit_cp_res.csv"))

#### Export Plots --------------------------------------------------------------
pdf(paste0(out_path,"any_visit_cp_plots.pdf"),onefile = TRUE)
for (i in 1:ceiling(nrow(final_results)/6)){
  tmp <- final_results %>%
    slice((i*6-5):(i*6)) %>%
    unnest(expected_trend) %>%
    inner_join(select(visit_counts,period,n),by = "period") %>%
    mutate(label = paste0("cp=",change_point," pb_cp=",cp_pred_bound," Method: ",method,"\n Periodicity: ",week_period,"; autoreg:",auto_reg))

  cps <- tmp %>%
    distinct(label,change_point,cp_pred_bound)

  p <- tmp %>%
    ggplot(aes(period,n)) +
    geom_point() +
    scale_x_reverse() +
    geom_line(aes(y = pred1),color = "red") +
    facet_wrap(~label, ncol = 2) +
    geom_vline(data = cps, aes(xintercept=change_point), linetype = 2) +
    geom_vline(data = cps, aes(xintercept=cp_pred_bound), linetype = 2, color = "green") +
    theme_bw()

  print(p)
}
dev.off()



##########################################
#### Find Change-point for ssd visits ####
##########################################

if (cond_name %in% codeBuildr::avail_ssd_codes(F)){
  ssd_codes <- bind_rows(codeBuildr::load_ssd_codes(cond_name) %>%
                           filter(type == "icd9") %>%
                           select(dx=code) %>%
                           mutate(icd_ver=9),
                         codeBuildr::load_ssd_codes(cond_name) %>%
                           filter(type == "icd10") %>%
                           select(dx=code) %>%
                           mutate(icd_ver=10))


  # load(paste0(in_path,"/all_dx_visits.RData"))

  #Change this to read from wherever the data is stored for you
  visit_counts <- all_dx_visits %>%
    inner_join(ssd_codes) %>%
    distinct(enrolid,days_since_index) %>%
    count(days_since_index) %>%
    filter(days_since_index<0) %>%
    mutate(period = -days_since_index) %>%
    select(period,n,days_since_index)


  max_days <- -1*min(visit_counts$days_since_index)

  ##################################################
  ### Calculate data for all changepoint methods ###
  ##################################################

  all_methods <- tibble(method = c("lm", "lm_quad", "lm_cube", "quad", "cube", "exp",
                                   # "spline",
                                   "pettitt", "cusum"
                                   # "MSE", "RMSE", "MAE", "MSLE", "RMSLE"
                                   )) %>%
    mutate(week_period = map(method, ~c(TRUE, FALSE))) %>% unnest(cols = week_period) %>%
    mutate(auto_reg = map(method, ~c(FALSE))) %>% unnest(cols = auto_reg) %>%
    # mutate(auto_reg = map(method, ~c(TRUE, FALSE))) %>% unnest(cols = auto_reg) %>%
    mutate(var = map(method, ~c("n"))) %>% unnest(cols = var)


  final_prior_visit_counts <- all_methods %>% mutate(counts = pmap(list(method),~filter(visit_counts,period<=max_days)))



  ###Fitting models###
  #run set change point function
  results <- final_prior_visit_counts %>%
    mutate(out=pmap(list(counts, method, week_period, auto_reg),
                    ~delaySim::find_change_point(data = ..1, method = ..2,
                                                 var_name = "n",
                                                 week_period = ..3, specify_cp = NULL,
                                                 auto_reg = ..4)))

  results <- results %>% mutate(expected_trend=map(out,~.$pred %>% select(period, t, pred1)))

  #extract change points
  final_results <- results %>%
    mutate(change_point=map(out, ~.$change_point$period)) %>%
    unnest(change_point)

  final_results <- final_results %>% select(-counts, -out)


  #################We now do a similar procedure, but for the prediction bound method#############################




  #Create a function to generate this data
  generate_pred_bound_results <- function(data, week_period, auto_reg , set_low){

    if(week_period | auto_reg){
      this_model <- "lm_period"
    } else{
      this_model <- "lm"
    }

    ## Prediction bound method with 90% bounds ##
    pred_boud_lm <- predictBoundCP::fit_cp_model(count_data = data %>%
                                                   mutate(days_since_dx = period * -1) %>%
                                                   select(days_since_dx, n),
                                                 low = set_low, high = max_days,
                                                 model = this_model,
                                                 level = 0.90,
                                                 plot = FALSE)

    #run set change point function to generate the results we seek
    results <- data %>%
      delaySim::find_change_point(method = "cusum",
                                  specify_cp = pred_boud_lm$cp,
                                  var_name = "n",
                                  week_period = week_period,
                                  auto_reg = auto_reg)

    #Change plots labels
    results$cp_plot$labels$title <- paste0("Method = pred_bound, week_effect = ",week_period,", auto_reg=",auto_reg)

    return(results)
  }

  ###Now, we prepare data and run this function
  # specify low to be 10% greater than cp found by cusum
  #run set pred bound function
  results_pred_bound <- final_prior_visit_counts %>% filter(method == "cusum") %>%
    mutate(low = pmap(list(week_period, auto_reg, var), ~round(final_results %>%
                                                                 filter(method == "cusum",
                                                                        week_period == ..1,
                                                                        auto_reg == ..2,
                                                                        var == ..3) %>% .$change_point * 1.10, 0)),
           out = pmap(list(counts, week_period, auto_reg, low),
                      ~generate_pred_bound_results(data = ..1,
                                                   week_period = ..2,
                                                   auto_reg = ..3,
                                                   set_low = ..4)))

  results <- results_pred_bound %>% mutate(expected_trend=map(out,~.$pred %>% select(period, t, pred1)))

  #extract change points
  final_results2 <- results %>%
    mutate(change_point=map(out, ~.$change_point$period)) %>%
    unnest(change_point)


  ###Aaron, I leave this in here if you want these plots, since they are
  ###already calculated by the find_change_point function automatically
  #Uncomment to include plots in the final_results object
  #create plots of change points
  #final_results2 <- final_results2 %>%
  #  mutate(graphs=map(out, ~.$cp_plot))

  final_results2 <- final_results2 %>% select(-out, -counts, -low) %>% mutate(method = "pred_bound")

  ########We append this data to our other datasets##########

  final_results <- rbind(final_results, final_results2)

  #Save these results
  second_pass_pred_bound <- function(pred_data,count_data){
    tmp <- pred_data %>%
      inner_join(count_data,by = "period") %>%
      arrange(period) %>%
      mutate(above_cp=cumsum(n>pred1)==row_number()) %>%
      filter(above_cp) %>%
      summarise(cp=max(period))

    as.integer(tmp$cp)
  }


  final_results <- final_results %>%
    mutate(cp_pred_bound = map_int(expected_trend,~second_pass_pred_bound(.,visit_counts)))

  #### Export Final change-points ------------------------------------------------
  final_results %>%
    select(method,periodicity=week_period,auto_reg,change_point,cp_pred_bound) %>%
    write_csv(paste0(out_path,"ssd_visit_cp_res.csv"))

  #### Export Plots --------------------------------------------------------------
  pdf(paste0(out_path,"ssd_visit_cp_plots.pdf"),onefile = TRUE)
  for (i in 1:ceiling(nrow(final_results)/6)){
    tmp <- final_results %>%
      slice((i*6-5):(i*6)) %>%
      unnest(expected_trend) %>%
      inner_join(select(visit_counts,period,n),by = "period") %>%
      mutate(label = paste0("cp=",change_point," pb_cp=",cp_pred_bound," Method: ",method,"\n Periodicity: ",week_period,"; autoreg:",auto_reg))

    cps <- tmp %>%
      distinct(label,change_point,cp_pred_bound)

    p <- tmp %>%
      ggplot(aes(period,n)) +
      geom_point() +
      scale_x_reverse() +
      geom_line(aes(y = pred1),color = "red") +
      facet_wrap(~label, ncol = 2) +
      geom_vline(data = cps, aes(xintercept=change_point), linetype = 2) +
      geom_vline(data = cps, aes(xintercept=cp_pred_bound), linetype = 2, color = "green") +
      theme_bw()

    print(p)
  }
  dev.off()

}

