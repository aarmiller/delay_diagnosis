
rm(list = ls())

########################
#### Initial Params ####
########################

true_cp <- -30



###############
#### Setup ####
###############

### Collect data to generate initial trends ------------------------------------

# Here we will start by collecting the basic counts of all visits before dengue
# diagnosis

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/dengue/dengue.db")

tm <- db %>% 
  tbl("tm") %>% 
  collect()

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()

visit_counts <- index_dates %>% 
  distinct(patient_id,index_date) %>% 
  inner_join(tm) %>% 
  mutate(days_since_index = svcdate-index_date) %>% 
  filter(between(days_since_index,-180,0)) %>% 
  filter(inpatient == 1 | ed == 1 | outpatient == 1 | obs_stay == 1) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index)

# add DOW factor
visit_counts <- visit_counts %>% 
  filter(days_since_index<0) %>% 
  mutate(dow = days_since_index %% 7) %>% 
  mutate(dow = as.factor(dow))

visit_counts %>% 
  filter(days_since_index<0) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  theme_bw() +
  geom_vline(aes(xintercept = -30), linetype = 2)


write_csv(visit_counts, file = "cp_approaches/simulate_trends/dengue_all_vis_counts.csv")

### Estimate basic trends ------------------------------------------------------

find_fits <- function(data,cp) {
  
  if(cp>0){ 
    stop("CP must be <0")
  } 
  
  # fit trends for linear, quadratic, and cubic
  expected_fit_lm <- lm(n~days_since_index + dow, data = filter(data,days_since_index < cp))
  expected_fit_quad <- lm(n~days_since_index+I(days_since_index^2) + dow, data = filter(data,days_since_index < cp))
  expected_fit_cube <- lm(n~days_since_index+I(days_since_index^2)+I(days_since_index^3) + dow, data = filter(data,days_since_index < cp))
  
  # add fitted values to data
  expected_fit_data <- data %>% 
    mutate(linear = predict(expected_fit_lm, newdata = .)) %>% 
    mutate(quad = predict(expected_fit_quad, newdata = .)) %>% 
    mutate(cube = predict(expected_fit_cube, newdata = .)) %>% 
    gather(key = model, value = pred, -days_since_index, -n, -dow)
  
  # Compute excess fits for data beyond CP
  excess_fits <- expected_fit_data  %>% 
    mutate(excess = n-pred) %>% 
    filter(days_since_index >= cp) %>%  
    group_by(model) %>% 
    nest() %>% 
    mutate(fit = map(data,~loess(excess~days_since_index, data = .))) %>% 
    mutate(fit = map(fit,~predict(., se = TRUE))) %>% 
    mutate(data = map2(data,fit, ~mutate(.x, excess_fit=.y$fit))) %>% 
    mutate(data = map2(data,fit, ~mutate(.x, excess_se=.y$se.fit))) %>% 
    select(model,data) %>% 
    unnest(data) %>% 
    ungroup() %>% 
    select(model,days_since_index,excess:excess_se)
  
  # generate empirical error terms for prior period
  prior_error_empirical <- expected_fit_data %>% 
    filter(days_since_index < cp) %>%
    mutate(error = pred-n) %>% 
    select(model,error) %>% 
    group_by(model) %>% 
    nest() %>% 
    deframe()
  
  prior_error_stats <- prior_error_empirical %>% 
    enframe(name = "model",value = "error") %>% 
    mutate(dist_params = map(error,~summarise(.,
                                              mean=mean(error),
                                              sd = sd(error)))) %>% 
    select(model,dist_params) %>% 
    unnest(dist_params)
  
  # final output
  out <- list(upper_bound = nrow(data),
              cp = cp,
              linear_model = coef(expected_fit_lm),
              quad_model = coef(expected_fit_quad),
              cube_model = coef(expected_fit_cube),
              expected_fits = expected_fit_data,
              excess_fits = excess_fits,
              prior_error_empirical = prior_error_empirical,
              prior_error_stats = prior_error_stats)
  
  return(out)
  
}


nrow(visit_counts)

fit_data <- find_fits(visit_counts,cp = -20)

# Plot estimated trends
fit_data$expected_fits %>% 
  ggplot(aes(x = days_since_index)) +
  geom_line(aes(y = pred, color = model), size = 1) +
  geom_point(aes(y = n)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_vline(aes(xintercept = fit_data$cp), linetype = 2)

# Plot excess fits (fitted by loess model)
fit_data$excess_fits %>% 
  mutate(pred_high = excess_fit + 1.96*excess_se,
         pred_low = excess_fit - 1.96*excess_se) %>% 
  ggplot(aes(days_since_index,excess_fit)) +
  geom_line() +
  geom_point(aes(y = excess)) +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha = 0.5) +
  facet_wrap(~model) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 1)


#####################################
#### Final Simulation Parameters ####
#####################################


tmp_fit_data <- find_fits(visit_counts,cp = -20)

tmp_fit_data$prior_error %>% 
  enframe(name = "model",value = "error") %>% 
  mutate(dist_params = map(error,~summarise(.,
                                            mean=mean(error),
                                            sd = sd(error)))) %>% 
  select(model,dist_params) %>% 
  unnest(dist_params)

tmp_fit_data$prior_error$linear %>% 
  ggplot(aes(error)) +
  geom_histogram()

tmp_fit_data$prior_error_stats %>% 
  mutate(error = map(sd,~rnorm(n = tmp_fit_data$upper_bound, mean = 0, sd =.))) %>% 
  unnest(error) %>% 
  select(model,error) %>% 
  mutate(days_since_index = rep(seq(from = -180, to = -1), times = 3))

?seq

tmp_fit_data$upper_bound

data <- tmp_fit_data

simulate_visits <- function(data, empirical_error = TRUE){
  
  ## Extract Expected ----------------------------------------------------------
  tmp_expected <- data$expected_fits %>% 
    select(model,days_since_index,dow,expected=pred)
  
  ## Draw Error ----------------------------------------------------------------
  if (empirical_error==TRUE){
    # Draw error term from empirical error values
    tmp_error <- map(data$prior_error_empirical, ~sample(.$error, size = data$upper_bound, replace = TRUE)) %>% 
      enframe(name = "model",value = "error") %>% 
      unnest("error") %>% 
      mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3)) 
  } else {
    # Draw error term from distribution
    tmp_error <- data$prior_error_stats %>% 
      mutate(error = map(sd,~rnorm(n = data$upper_bound, mean = 0, sd =.))) %>% 
      unnest(error) %>% 
      select(model,error) %>% 
      mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3))
  }
  
  # Draw Excess visits during delay window
  tmp_excess <- data$excess_fits %>% 
    mutate(excess = map2_dbl(excess_fit,excess_se,~rnorm(1,mean = .x, sd = .y))) %>% 
    select(model,days_since_index, excess)
  
  out <- tmp_expected %>% 
    inner_join(tmp_error, by = join_by(model, days_since_index)) %>% 
    left_join(tmp_excess, by = join_by(model, days_since_index)) %>% 
    mutate(excess = replace_na(excess, 0L)) %>% 
    mutate(simulated = expected+excess+error) %>% 
    select(model,days_since_index,dow,simulated)
  
  return(out)
}


simulate_visits(tmp_fit_data,FALSE)

run_sim <- function(data, n_trials = 100, empirical_error = TRUE){
  tibble(trial = 1:n_trials) %>% 
    mutate(simulated_data = map(trial,~simulate_visits(data = data, empirical_error = empirical_error)))
}

#############################
#### Evaluate Simulation ####
#############################



## Compare using Empirical versus distributional error terms -------------------
tmp <- bind_rows(run_sim(data = tmp_fit_data) %>% 
            mutate(error_type = "empirical"),
          run_sim(data = tmp_fit_data, empirical_error = FALSE) %>% 
            mutate(error_type = "distributional"))

tmp %>% 
  unnest(simulated_data) %>% 
  ggplot(aes(days_since_index,simulated, group = trial, color=error_type)) +
  geom_line() +
  facet_wrap(~model)

tmp %>% 
  unnest(simulated_data) %>% 
  ggplot(aes(days_since_index,simulated, group = trial)) +
  geom_line() +
  facet_grid(model~error_type)


## Compare using different cutoffs ---------------------------------------------

n_trials <- 500

tmp <- tibble(cp = c(-10,-20,-30,-40)) %>% 
  mutate(sim_data = map(cp,~find_fits(visit_counts, cp = .) %>% 
                          run_sim(n_trials = n_trials))) %>% 
  unnest(sim_data)


# Compute MSE over entire
tmp %>% 
  unnest(simulated_data) %>% 
  inner_join(select(visit_counts, days_since_index,n), by = "days_since_index") %>% 
  group_by(cp,model) %>% 
  summarise(mse_overall = mean((n-simulated)^2))

# Compute MSE over respective prior
tmp %>% 
  unnest(simulated_data) %>% 
  inner_join(select(visit_counts, days_since_index,n), by = "days_since_index") %>% 
  filter(days_since_index<cp) %>% 
  group_by(cp,model) %>% 
  summarise(mse_prior = mean((n-simulated)^2)) %>% 
  ungroup() %>% 
  spread(key = model, value = mse_prior)

# Compute MSE over similar period
tmp %>% 
  unnest(simulated_data) %>% 
  inner_join(select(visit_counts, days_since_index,n), by = "days_since_index") %>% 
  filter(days_since_index< -40) %>% 
  group_by(cp,model) %>% 
  summarise(mse_prior = mean((n-simulated)^2)) %>% 
  ungroup() %>% 
  spread(key = model, value = mse_prior)


# Linear Models
bind_rows(tmp %>% 
            unnest(simulated_data) %>% 
            filter(model == "linear"),
          visit_counts %>% 
            mutate(cp = map(days_since_index, ~c(-10,-20,-30,-40))) %>% 
            unnest(cp)) %>% 
  ggplot(aes(days_since_index,simulated, group = trial)) +
  geom_line() +
  geom_point(aes(y=n), color = "red") +
  geom_vline(aes(xintercept = cp), linetype = 2) +
  facet_wrap(~cp)


# Quadratic Models
bind_rows(tmp %>% 
            unnest(simulated_data) %>% 
            filter(model == "quad"),
          visit_counts %>% 
            mutate(cp = map(days_since_index, ~c(-10,-20,-30,-40))) %>% 
            unnest(cp)) %>% 
  ggplot(aes(days_since_index,simulated, group = trial)) +
  geom_line() +
  geom_point(aes(y=n), color = "red") +
  geom_vline(aes(xintercept = cp), linetype = 2) +
  facet_wrap(~cp)

# Cubic Models
bind_rows(tmp %>% 
            unnest(simulated_data) %>% 
            filter(model == "cube"),
          visit_counts %>% 
            mutate(cp = map(days_since_index, ~c(-10,-20,-30,-40))) %>% 
            unnest(cp)) %>% 
  ggplot(aes(days_since_index,simulated, group = trial)) +
  geom_line() +
  geom_point(aes(y=n), color = "red") +
  geom_vline(aes(xintercept = cp), linetype = 2) +
  facet_wrap(~cp)


##############
#### Idea ####
##############

# Empirical Specification test

n_trials <- 500

tmp <- tibble(cp = -8:-40) %>% 
  mutate(fits = map(cp,~find_fits(visit_counts,cp = .))) %>% 
  mutate(sim_data = map(fits,~run_sim(data = ., n_trials = n_trials, empirical_error = FALSE)))


spec_eval <- tmp %>% 
  select(cp,sim_data) %>% 
  unnest(sim_data) %>% 
  unnest(simulated_data) %>% 
  group_by(cp,model,days_since_index) %>% 
  summarise(sim_high = quantile(simulated, probs = c(0.99)),
            sim_low = quantile(simulated, probs = c(0.01))) %>% 
  ungroup() %>% 
  inner_join(select(visit_counts,days_since_index,n), by = "days_since_index") %>% 
  group_by(cp,model) %>% 
  summarise(pct_contained=sum(n<sim_high & n>sim_low)/n()) %>% 
  arrange(desc(pct_contained)) %>% 
  ungroup()

find_fits(visit_counts,cp = -18)$expected_fits %>% 
  filter(model == "cube") %>% 
  ggplot(aes(days_since_index)) +
  geom_point(aes(y = n)) +
  geom_line(aes(y = pred)) +
  geom_vline(aes(xintercept = -18))

find_fits(visit_counts,cp = -17)$expected_fits %>% 
  filter(model == "quad") %>% 
  ggplot(aes(days_since_index)) +
  geom_point(aes(y = n)) +
  geom_line(aes(y = pred)) +
  geom_vline(aes(xintercept = -17))

tmp %>% 
  mutate(expected_fits = map(fits, ~.$expected_fits)) %>% 
  select(cp,expected_fits) %>% 
  unnest(expected_fits) 

# Compare Fit of Selected Models
spec_eval %>% 
  filter(pct_contained>.99) %>% 
  inner_join(tmp %>% 
               mutate(expected_fits = map(fits, ~.$expected_fits)) %>% 
               select(cp,expected_fits) %>% 
               unnest(expected_fits)) %>% 
  mutate(label = paste0("Model: ", model,", CP: ",cp)) %>% 
  ggplot(aes(days_since_index)) +
  geom_point(aes(y = n)) +
  geom_line(aes(y = pred), color = "red", size = 1) +
  facet_wrap(~label) +
  geom_vline(aes(xintercept = cp), linetype = 2)

# Compare Implied Excess

spec_eval %>% 
  filter(pct_contained>.99) %>% 
  inner_join(tmp %>% 
               mutate(excess_fits = map(fits, ~.$excess_fits)) %>% 
               select(cp,excess_fits) %>% 
               unnest(excess_fits)) %>% 
  mutate(label = paste0("Model: ", model,", CP: ",cp)) %>% 
  mutate(pred_high = excess_fit + 1.96*excess_se,
         pred_low = excess_fit - 1.96*excess_se) %>% 
  ggplot(aes(days_since_index)) +
  geom_point(aes(y = excess)) +
  facet_wrap(~label) +
  geom_line(aes(y = excess_fit)) +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha = 0.5) +
  geom_hline(aes(yintercept = 0))



##############################
#### Add Bootstrapp Layer ####
##############################

tmp_tm <- index_dates %>% 
  distinct(patient_id,index_date) %>% 
  inner_join(tm) %>% 
  mutate(days_since_index = svcdate-index_date) %>% 
  filter(between(days_since_index,-180,-1)) %>% 
  filter(inpatient == 1 | ed == 1 | outpatient == 1 | obs_stay == 1) %>% 
  distinct(patient_id,days_since_index)

ids <- distinct(index_dates,patient_id,index_date)

draw_boot_counts <- function(){
  sample_n(ids, n(),replace = TRUE) %>% 
    mutate(boot_id = row_number()) %>% 
    inner_join(tmp_tm,relationship = "many-to-many",by = join_by(patient_id)) %>% 
    distinct(boot_id,patient_id,days_since_index) %>% 
    count(days_since_index) %>% 
    mutate(dow = days_since_index %% 7) %>% 
    mutate(dow = as.factor(dow))
}

boot_counts <- tibble(bootstrap = 1:10) %>% 
  mutate(visit_counts = map(bootstrap,~draw_boot_counts()))

n_trials <- 500

sim_res <- list()

for (i in 1:nrow(boot_counts)){
  
  print(paste0("Bootstrap # ", i))
  
  tmp_visit_counts <- boot_counts$visit_counts[[i]]
  
  tmp <- tibble(cp = -8:-40) %>% 
    mutate(fits = map(cp,~find_fits(tmp_visit_counts,cp = .))) %>% 
    mutate(sim_data = map(fits,~run_sim(data = ., n_trials = n_trials, empirical_error = FALSE)))
  
  
  sim_res[[i]] <- tmp
  

}

sim_res

boot_counts$visit_counts[[1]]
sim_res[[1]]

evaluate_specification <- function(visit_count_data,sim_res_data){
  sim_res_data %>% 
    select(cp,sim_data) %>% 
    unnest(sim_data) %>% 
    unnest(simulated_data) %>% 
    group_by(cp,model,days_since_index) %>% 
    summarise(sim_high = quantile(simulated, probs = c(0.99)),
              sim_low = quantile(simulated, probs = c(0.01))) %>% 
    ungroup() %>% 
    inner_join(select(visit_count_data,days_since_index,n), by = "days_since_index") %>% 
    group_by(cp,model) %>% 
    summarise(pct_contained=sum(n<sim_high & n>sim_low)/n()) %>% 
    arrange(desc(pct_contained)) %>% 
    ungroup()
}

evaluate_specification(boot_counts$visit_counts[[1]],sim_res[[1]])

boot_spec_eval <- map2(boot_counts$visit_counts,sim_res,evaluate_specification)

boot_spec_results <- bind_rows(boot_spec_eval) %>% 
  group_by(cp,model) %>% 
  summarise(pct_contained = median(pct_contained)) %>% 
  ungroup() %>% 
  arrange(desc(pct_contained))

boot_excess_fits <- map(sim_res, ~select(.,cp,fits)) %>% 
  enframe(name = "bootstrap") %>% 
  unnest(value) %>% 
  mutate(excess_fits = map(fits, ~.$excess_fits)) %>% 
  select(bootstrap,cp,excess_fits) %>% 
  unnest(excess_fits)

boot_spec_results

boot_excess_fits %>% 
  filter(cp == -20,
         model == "linear") %>% 
  mutate(pred_high = excess_fit + 1.96*excess_se,
         pred_low = excess_fit - 1.96*excess_se) %>% 
  ggplot(aes(days_since_index)) +
  geom_point(aes(y = excess)) +
  facet_wrap(~bootstrap) +
  geom_line(aes(y = excess_fit)) +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_high), alpha = 0.5) +
  geom_hline(aes(yintercept = 0))




34




#########################
#### Simulate Trends ####
#########################

base_data <- fit_data %>% 
  inner_join(params$dow_effects) %>% 
  select(-n)

### Functions ------------------------------------------------------------------

draw_excess <- function(data){
  mutate(data,excess = map2_dbl(excess_fit,excess_se,~rnorm(1,mean = .x, sd = .y))) %>% 
    select(days_since_index, excess)}

# make_expected_data <- function(dow_shift = 0) {
#   select(fit_data,days_since_index) %>% 
#     mutate(dow = as.factor((days_since_index+dow_shift) %% 7)) %>% 
#     inner_join(dow_effects,by = join_by(dow)) %>% 
#     mutate(error = rnorm(n(), mean = 0, sd = error_sigma)) %>% 
#     mutate(expected = intercept + effect + days_since_index*trend + error) %>% 
#     select(days_since_index,dow,expected)
# }



simulate_visits <- function(params){
  
  error_terms <- map(params$prior_error,~sample(.$error, size = 180, replace = TRUE)) %>% 
    enframe(name = "model",value = "error") %>% 
    unnest("error") %>% 
    mutate(days_since_index = rep((-180):(-1),time = 3))
  
  tmp1 <- base_data %>% 
    mutate(linear = params$linear_terms[1] + params$linear_terms[2]*days_since_index + dow_effect) %>% 
    mutate(quad = params$quad_terms[1] + params$quad_terms[2]*days_since_index + params$quad_terms[3]*(days_since_index)^2 + dow_effect) %>% 
    mutate(cube = params$cube_terms[1] + params$cube_terms[2]*days_since_index + params$cube_terms[3]*(days_since_index)^2 + 
             params$cube_terms[4]*(days_since_index)^3 + dow_effect) %>% 
    gather(key = model, value = expected, -days_since_index, -dow, -dow_effect) %>% 
    left_join(error_terms, by = c("model", "days_since_index"))
  
  
  # draw excess visits during delay window
  tmp_excess <- rename(draw_excess(params$excess_fits[[1]]),linear=excess) %>% 
    inner_join(rename(draw_excess(params$excess_fits[[2]]),quad=excess), by = join_by(days_since_index)) %>% 
    inner_join(rename(draw_excess(params$excess_fits[[3]]),cube=excess), by = join_by(days_since_index)) %>% 
    gather(key = model, value = excess, -days_since_index)
  
  out <- tmp1 %>% 
    left_join(tmp_excess,by = join_by(days_since_index, model)) %>% 
    mutate(excess = replace_na(excess, 0L)) %>% 
    mutate(observed = expected+excess+error) %>% 
    select(days_since_index,dow,model,observed)
  
  return(out)
}


simulate_visits() %>% 
  ggplot(aes(days_since_index,observed, color = model)) +
  geom_line()


#######################
#### Simulate Data ####
#######################

n_trials <- 1000

simulated_data <- tibble(trial = 1:n_trials) %>% 
  mutate(data = map(trial,~simulate_visits()))

bind_rows(simulated_data %>% 
            unnest(data),
          mutate(fit_data,model = map(days_since_index, ~c("linear","quad","cube"))) %>% 
                 unnest(model)) %>% 
  ggplot(aes(days_since_index, observed, group = trial)) +
  geom_line( alpha = 0.1) +
  facet_wrap(~model) +
  geom_point(aes(y = n), color = "red") +
  theme_bw()



#############
#### OLD ####
#############


base_data <- fit_data %>% 
  rename(true_n = n) %>% 
  inner_join(params$dow_effects)  





simulate_visits <- function(){
  
  tmp1 <- base_data %>% 
    mutate(error = rnorm(n(), mean = 0, sd = params$prior_se)) %>% 
    mutate(linear_expected = params$linear_terms[1] + params$linear_terms[2]*days_since_index + dow_effect + error) %>% 
    mutate(quad_expected = params$quad_terms[1] + params$quad_terms[2]*days_since_index + params$quad_terms[3]*(days_since_index)^2 + dow_effect + error) %>% 
    mutate(cube_expected = params$cube_terms[1] + params$cube_terms[2]*days_since_index + params$cube_terms[3]*(days_since_index)^2 + 
             params$cube_terms[4]*(days_since_index)^3 + dow_effect + error)
  
  
  # draw excess visits during delay window
  tmp_excess <- rename(draw_excess(params$excess_fits[[1]]),linear_excess=excess) %>% 
    inner_join(rename(draw_excess(params$excess_fits[[2]]),quad_excess=excess), by = join_by(days_since_index)) %>% 
    inner_join(rename(draw_excess(params$excess_fits[[3]]),cube_excess=excess), by = join_by(days_since_index))
  
  out <- tmp1 %>% 
    left_join(tmp_excess,by = join_by(days_since_index)) %>% 
    mutate_at(vars(linear_excess:cube_excess), ~replace_na(.,0)) %>% 
    mutate(linear_obs = linear_expected + linear_excess,
           quad_obs = quad_expected + quad_excess,
           cube_obs = cube_expected + cube_excess) %>% 
    # select(days_since_index,true_n,linear_obs:cube_obs) 
    select(days_since_index,dow,linear_obs:cube_obs)
  
  return(out)
}

simulate_visits() %>% 
  gather(key = key, value = value, -days_since_index, -dow) %>% 
  ggplot(aes(x = days_since_index)) +
  geom_line(aes(y = value, color = key))


#######################################
#### Evaluate Chane-point approach ####
#######################################


simulated_data <- tibble(trial = 1:1000) %>% 
  mutate(data = map(trial,~simulate_visits()))

bind_rows(simulated_data %>% 
              unnest(data),base_data) %>% 
  ggplot(aes(days_since_index, cube_obs, group = trial)) +
  geom_line(alpha = .2) +
  geom_point(aes(x=days_since_index,y=true_n), color = "red") +
  theme_bw()


find_implied_cp <- function(data){
  tmp <- data %>% 
    arrange(desc(days_since_index)) %>% 
    mutate(excess_obs = cumsum(linear_obs>pred)) %>% 
    filter(excess_obs == row_number()) 
  
  
  min(tmp$days_since_index)
    
}



simulated_data$data[[1]] %>% 
  mutate(n = linear_obs) %>% 
  fit_model(cp = 40, return_fit = TRUE)


compute_misses <- function(data, cp){
  data %>% 
    filter(days_since_index>cp) %>% 
    summarise(misses = sum(n-pred)) %>% 
    .$misses
}

tmp <- simulated_data %>% 
  mutate(data = map(data,~mutate(.,n = linear_obs))) %>% 
  mutate(data = map(data,~fit_prior_period(., cp = -40))) %>% 
  mutate(implied_cp = map_int(data,find_implied_cp)) %>% 
  mutate(num_miss = map2_dbl(data,implied_cp,compute_misses))

tmp


# now compute for the overall sample

base_data %>% 
  mutate(n = true_n) %>% 
  fit_prior_period(cp = -30) %>% 
  filter(days_since_index > -30) %>% 
  mutate(misses = n-pred) %>% 
  filter(misses>0) %>% 
  summarise(misses = sum(misses))

3469




tmp %>% 
  select(trial,implied_cp,num_miss) %>% 
  ggplot(aes(implied_cp)) +
  geom_histogram()


tmp %>% 
  select(trial,implied_cp,num_miss) %>% 
  ggplot(aes(num_miss)) +
  geom_histogram() +
  geom_vline(aes(xintercept = 3469), color = "red") +
  geom_vline(aes(xintercept = 3466), color = "black") +
  xlab("Number of Missed Visits") +
  theme_bw()

tmp %>% 
  select(trial,implied_cp,num_miss) %>% 
  summarise(mean(num_miss))


tmp$data[[1]] %>% 
  filter(days_since_index>-23) %>% 
  summarise(misses = sum(n-pred)) %>% 
  .$misses


tmp$data[[2]] %>% 
  arrange(desc(days_since_index)) %>% 
  mutate(excess_obs = cumsum(linear_obs>pred)) %>% 
  filter(excess_obs)

tmp$data[[2]] %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  geom_line(aes(y = pred), color = "red")


fit_prior_period <- function(data,cp,model = "linear"){
  tmp <- data %>% 
    filter(days_since_index <= cp)
  
  if (model=="linear"){
    fit <- lm(n~days_since_index+dow, data = tmp)
  } else if (model=="quad"){
    fit <- lm(n~poly(days_since_index,2)+dow, data = tmp)
  } else if (model=="cube") {
    fit <- lm(n~poly(days_since_index,3)+dow, data = tmp)
  }
  
  out <- data %>% 
    mutate(pred = predict(fit, newdata=data))
  
  return(out)
  
}


tmp <- simulated_data %>% 
  mutate(data = map(data,~mutate(.,n = linear_obs))) %>% 
  mutate(data = map(data,~fit_prior_period(data = ., cp = -40)))



tmp$data[[2]] %>% 
  arrange(desc(days_since_index)) %>% 
  mutate(excess_obs = cumsum(linear_obs>pred)) %>% 
  filter(excess_obs == row_number()) %>% 
  summarise(implied_cp = min(days_since_index))



tmp$data[[2]] %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  geom_line(aes(y = pred), color = "red")




fit_model <- function(data,cp,model="linear",return_fit=FALSE){
  
  tmp <- data %>% 
    mutate(before_cp = days_since_index <= cp)
  
  if (model=="linear"){
    fit <- lm(n~days_since_index*before_cp+dow, data = tmp)
  } else if (model=="quad"){
    fit <- lm(n~poly(days_since_index,2)*before_cp+dow, data = tmp)
  } else if (model=="cube") {
    fit <- lm(n~poly(days_since_index,3)*before_cp+dow, data = tmp)
  }
  
  rmse <- sqrt(mean(fit$residuals^2))
  
  if (return_fit){
    return(list(fit=fit,
                rmse=rmse,
                preds = predict(fit, new_data = data)))
  } else {
    return(rmse)
  }
  
}

find_cp <- function(data,cp_range,model="linear"){
  out1 <- tibble(cp = cp_range) %>% 
    mutate(rmse = map_dbl(cp,~fit_model(data = data,
                                        cp = .,
                                        model = model))) %>% 
    filter(rmse==min(rmse))
  
  out2 <- fit_model(data = data,
                    cp = out1$cp,
                    model = model,
                    return_fit = TRUE)
  
  out3 <- data %>% 
    mutate(before_cp = TRUE) %>%
    mutate(pred1=predict(out2$fit,newdata=.)) %>%
    mutate(before_cp = days_since_index<=out1$cp) %>% 
    mutate(pred2=predict(out2$fit,newdata=.)) %>% 
    select(days_since_index,n,pred1,pred2)
  
  return(list(cp=out1$cp,
              fit=out2$fit,
              pred=out3))
}

rename(simulated_data$data[[4]], n = linear_obs) %>% 
  fit_model(cp = -30, model = "linear",TRUE)


rename(simulated_data$data[[4]], n = linear_obs) %>% 
  find_cp(data = ., cp_range = -50:(-5), model = "linear")

tmp <- simulated_data %>% 
  mutate(data = map(data, ~mutate(., n = linear_obs))) %>% 
  mutate(fit = map(data,~find_cp(data=.,cp_range = -50:(-5), model = "linear")))

tmp

tmp <- tmp %>% 
  mutate(found_cp = map_int(fit,~.$cp)) %>% 
  mutate(found_slope = map_dbl(fit,~coef(.$fit)["days_since_index:before_cpTRUE"])) 

tmp %>% 
  mutate(found_cp = map_int(fit,~.$cp)) %>% 
  select()

fit_prior <- function(data,cp, model = "linear"){
  tmp <- data %>% 
    filter(days_since_index <= cp)
  
  if (model=="linear"){
    fit <- lm(n~days_since_index+dow, data = tmp)
  } else if (model=="quad"){
    fit <- lm(n~poly(days_since_index,2)+dow, data = tmp)
  } else if (model=="cube") {
    fit <- lm(n~poly(days_since_index,3)+dow, data = tmp)
  }
}

tmp2 <- tmp %>% 
  select(trial,data,found_cp) %>% 
  mutate(fit = map2(data,found_cp, ~fit_prior(.x,.y))) %>% 
  mutate(slope = map_dbl(fit,~coef(.)["days_since_index"]))

tmp2$fit[1]

fit_shifted_cp <- function(shift){
  tmp %>% 
    select(trial,data, found_cp) %>% 
    mutate(fit = map2(data,found_cp, ~fit_prior(.x,.y+shift))) %>% 
    mutate(slope = map_dbl(fit,~coef(.)["days_since_index"])) %>% 
    select(trial,slope) 
}

fit_shifted_cp(-5)


shift_res <- tibble(shift = -30:-1) %>% 
  mutate(res = map(shift,fit_shifted_cp))


shift_res <- bind_rows(tmp %>% select(trial,data,found_cp) %>% 
            mutate(fit = map2(data,found_cp, ~fit_prior(.x,.y))) %>% 
            mutate(slope = map_dbl(fit,~coef(.)["days_since_index"])) %>% 
            select(trial, slope) %>% 
            mutate(shift = 0),
          shift_res %>% unnest(res))

shift_centers <- shift_res %>% 
  mutate(shift = -shift) %>% 
  group_by(shift) %>% 
  summarise(mean_slope = mean(slope),
            median_slope = median(slope))

shift_res %>% 
  mutate(shift = -shift) %>% 
  ggplot(aes(slope)) +
  geom_histogram() +
  geom_vline(aes(xintercept = params$linear_terms[2]), color = "red") +
  geom_vline(aes(xintercept = mean_slope), data = shift_centers) +
  facet_wrap(~shift)

tmp2$fit[[1]]

tmp %>% 
  mutate(found_cp = map_int(fit,~.$cp)) %>% 
  mutate(found_slope = map_dbl(fit,~coef(.$fit)["days_since_index:before_cpTRUE"])) %>% 
  select(trial,data,found_cp,found_slope) %>% 
  mutate(fit_5 = map2(data,found_cp,~fit_model(data = .x, cp = .y+5, model = "linear", TRUE)))


coef(tmp$fit[[1]]$fit)["period"]

params$linear_terms

tmp



excess_fits <- filter(fit_data, days_since_index > -30) %>% 
  select(days_since_index) %>% 
  mutate(excess_fit = predict(miss_fit, se = TRUE)$fit,
         excess_se = predict(miss_fit, se = TRUE)$se.fit)



expected_fit_lm <- lm(n~days_since_index + dow, data = filter(fit_data,days_since_index < -30))
expected_fit_quad <- lm(n~days_since_index+I(days_since_index^2) + dow, data = filter(fit_data,days_since_index < -30))
expected_fit_cube <- lm(n~days_since_index+I(days_since_index^2)+I(days_since_index^3) + dow, data = filter(fit_data,days_since_index < -30))
expected_fit_exp <- lm(n~exp(days_since_index) + dow, data = filter(fit_data,days_since_index < -30))

coef(expected_fit_lm)
coef(expected_fit_quad)
coef(expected_fit_cube)
100*coef(expected_fit_cube)[4]
coef(expected_fit_exp)[1]

tmp <- fit_data %>% 
  mutate(pred = predict(expected_fit_quad, newdata = .)) 

tmp$pred[1]

fit_data <- fit_data %>% 
  mutate(pred1 = predict(expected_fit,newdata = fit_data)) %>% 
  mutate(excess=n-pred1)


fit_data %>%
  ggplot(aes(days_since_index, n)) +
  geom_point() +
  geom_line(aes(y = pred1))


fit_data %>%
  ggplot(aes(days_since_index, excess)) +
  geom_line()

miss_fit <- lm(excess~poly(days_since_index,4), data = filter(fit_data, days_since_index > -30))
miss_fit <- loess(excess ~ days_since_index,data = filter(fit_data, days_since_index > -30))


filter(fit_data, days_since_index > -30) %>% 
  mutate(pred_excess = predict(miss_fit)) %>% 
  ggplot(aes(days_since_index, excess)) +
  geom_point() +
  geom_line(aes(y = pred_excess))

predict(miss_fit, se = TRUE)$fit
predict(miss_fit, se = TRUE)$se.fit

excess_fits <- filter(fit_data, days_since_index > -30) %>% 
  select(days_since_index) %>% 
  mutate(excess_fit = predict(miss_fit, se = TRUE)$fit,
         excess_se = predict(miss_fit, se = TRUE)$se.fit)

draw_excess(params$excess_fits[[2]]) 





draw_excess <- function(data){ mutate(data,excess = map2_dbl(excess_fit,excess_se,~rnorm(1,mean = .x, sd = .y))) }



tibble(trial = 1:100) %>% 
  mutate(res = map(trial,~draw_excess())) %>% 
  unnest(res) %>% 
  ggplot(aes(days_since_index,n,group = trial)) +
  geom_line()


true_cp <- -30

intercept <- 270

dow_effects <- tibble(dow = as.factor(0:6),
                      effect = c(0, -52, -63, -67, -74, -71, -45))

trend <- 0.15

error_sigma <- 15

make_expected_data <- function(dow_shift = 0) {
  select(fit_data,days_since_index) %>% 
    mutate(dow = as.factor((days_since_index+dow_shift) %% 7)) %>% 
    inner_join(dow_effects,by = join_by(dow)) %>% 
    mutate(error = rnorm(n(), mean = 0, sd = error_sigma)) %>% 
    mutate(expected = intercept + effect + days_since_index*trend + error) %>% 
    select(days_since_index,dow,expected)
}


make_expected_data(4) %>% 
  left_join(draw_excess()) %>% 
  mutate(excess = replace_na(excess,0)) %>% 
  mutate(observed = expected + excess) %>% 
  ggplot(aes(days_since_index,observed)) +
  geom_line()


select(fit_data,days_since_index) %>% 
  mutate(dow = (days_since_index+5) %% 7) %>% 
  filter(between(days_since_index,-180,-25))


tmp <- make_expected_data(4) %>% 
  left_join(draw_excess(),by = join_by(days_since_index)) %>% 
  mutate(excess = replace_na(excess,0)) %>% 
  mutate(observed = expected + excess) %>% 
  filter(days_since_index<=-25)


coef(lm(observed~days_since_index+dow, data = tmp))["days_since_index"]

find_fit <- function(cutoff,dow_shift = 0){
  tmp <- make_expected_data(dow_shift) %>% 
    left_join(draw_excess(),by = join_by(days_since_index)) %>% 
    mutate(excess = replace_na(excess,0)) %>% 
    mutate(observed = expected + excess) %>% 
    filter(days_since_index <= cutoff)
  
  
  coef(lm(observed~days_since_index+dow, data = tmp))["days_since_index"]
}

find_fit(-30)


tibble(trial = 1:1000) %>% 
  mutate(est = map_dbl(trial,~find_fit(-25))) %>% 
  ggplot(aes(est)) +
  geom_histogram() +
  geom_vline(aes(xintercept = trend))

sim_res <- tibble(dow_shift = 0:6) %>% 
  mutate(trial = map(dow_shift,~1:1000)) %>% 
  unnest(trial) %>% 
  mutate(est = map_dbl(dow_shift,~find_fit(-25, dow_shift = .)))


sim_res %>% 
  ggplot(aes(est)) +
  geom_histogram() +
  geom_vline(aes(xintercept = trend)) +
  facet_wrap(~dow_shift)


sim_res %>% 
  group_by(dow_shift) %>% 
  summarise(median(est))


sim_res <- tibble(cutoff = -40:-20) %>% 
  mutate(trial = map(cutoff,~1:1000)) %>% 
  unnest(trial) %>% 
  mutate(est = map_dbl(cutoff,~find_fit(cutoff = .)))

sim_res %>% 
  group_by(cutoff) %>% 
  mutate(median_est = median(est)) %>% 
  ggplot(aes(est)) +
  geom_histogram() +
  facet_wrap(~cutoff) +
  geom_vline(aes(xintercept = trend), color = "black") +
  geom_vline(aes())
  

         