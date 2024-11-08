



rm(list = ls())


visit_counts <- read_csv("cp_approaches/simulate_trends/dengue_all_vis_counts.csv") %>% 
  mutate(dow = as.factor(dow))

### Fit Trends -----------------------------------------------------------------


find_fits <- function(data,cp) {
  
  if(cp>0){ 
    stop("CP must be <0")
  } 
  
  # fit trends for linear, quadratic, and cubic models -------------------------
  expected_fit_lm <- lm(n~days_since_index + dow, data = filter(data,days_since_index < cp))
  expected_fit_quad <- lm(n~days_since_index+I(days_since_index^2) + dow, data = filter(data,days_since_index < cp))
  expected_fit_cube <- lm(n~days_since_index+I(days_since_index^2)+I(days_since_index^3) + dow, data = filter(data,days_since_index < cp))
  
  # add fitted values to data along with SE for prediction bound ---------------
  expected_fit_data <- data %>% 
    mutate(linear = predict(expected_fit_lm, newdata = .)) %>% 
    mutate(linear_se = predict(expected_fit_lm, newdata = ., se =TRUE)$se.fit) %>% 
    mutate(quad = predict(expected_fit_quad, newdata = .)) %>% 
    mutate(quad_se = predict(expected_fit_quad, newdata = ., se =TRUE)$se.fit) %>% 
    mutate(cube = predict(expected_fit_cube, newdata = .)) %>% 
    mutate(cube_se = predict(expected_fit_cube, newdata = ., se =TRUE)$se.fit) %>% 
    gather(key = model, value = pred, -days_since_index, -n, -dow) %>% 
    mutate(pred_type = ifelse(str_detect(model,"se"),"pred_se","pred")) %>% 
    mutate(model = str_remove(model,"_se")) %>% 
    spread(key = pred_type, value = pred) %>% 
    arrange(model,days_since_index)
  
  # Compute excess fits for data beyond CP -------------------------------------
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
  
  # generate empirical error terms for prior period ----------------------------
  # prior_error_empirical <- expected_fit_data %>% 
  #   filter(days_since_index < cp) %>%
  #   mutate(error = pred-n) %>% 
  #   select(model,error) %>% 
  #   group_by(model) %>% 
  #   nest() %>% 
  #   deframe()
  
  # prior_error_stats <- prior_error_empirical %>% 
  #   enframe(name = "model",value = "error") %>% 
  #   mutate(dist_params = map(error,~summarise(.,
  #                                             mean=mean(error),
  #                                             sd = sd(error)))) %>% 
  #   select(model,dist_params) %>% 
  #   unnest(dist_params)
  
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


### Simulate data --------------------------------------------------------------

# Three types of error

simulate_visits <- function(data, error_type = "pred"){
  
  ## Extract Expected ----------------------------------------------------------
  tmp_expected <- data$expected_fits %>% 
    select(model,days_since_index,dow,expected=pred, expected_se=pred_se) %>% 
    left_join(select(data$excess_fits,-excess),by = join_by(model, days_since_index)) %>% 
    mutate_at(vars(excess_fit:excess_se), ~replace_na(.,0)) %>% 
    mutate(draw_mean = expected+excess_fit,
           draw_se = expected_se+excess_se)
  
  # ## Draw Error ----------------------------------------------------------------
  # if (empirical_error==TRUE){
  #   # Draw error term from empirical error values
  #   tmp_error <- map(data$prior_error_empirical, ~sample(.$error, size = data$upper_bound, replace = TRUE)) %>% 
  #     enframe(name = "model",value = "error") %>% 
  #     unnest("error") %>% 
  #     mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3)) 
  # } else {
  #   # Draw error term from distribution
  #   tmp_error <- data$prior_error_stats %>% 
  #     mutate(error = map(sd,~rnorm(n = data$upper_bound, mean = 0, sd =.))) %>% 
  #     unnest(error) %>% 
  #     select(model,error) %>% 
  #     mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3))
  # }
  
  # # Draw Excess visits during delay window
  # tmp_excess <- data$excess_fits %>% 
  #   mutate(excess = map2_dbl(excess_fit,excess_se,~rnorm(1,mean = .x, sd = .y))) %>% 
  #   select(model,days_since_index, excess)
  # 
  # out <- tmp_expected %>% 
  #   inner_join(tmp_error, by = join_by(model, days_since_index)) %>% 
  #   left_join(tmp_excess, by = join_by(model, days_since_index)) %>% 
  #   mutate(excess = replace_na(excess, 0L)) %>% 
  #   mutate(simulated = expected+excess+error) %>% 
  #   select(model,days_since_index,dow,simulated)
  
  # Draw visits
  out <- tmp_expected %>% 
    mutate(simulated = map2_dbl(draw_mean,draw_se,~rnorm(1,mean = .x, sd = .y))) %>% 
    select(model,days_since_index,dow,simulated)
  
  return(out)
}








find_fits(visit_counts,-30)


data <- find_fits(visit_counts,-20)

find_fits(visit_counts,-17) %>% 
  simulate_visits() %>% 
  ggplot(aes(days_since_index,simulated, color = model)) +
  geom_line()

visit_counts %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line()


simulate_visits <- function(data){
  
  ## Extract Expected ----------------------------------------------------------
  tmp_expected <- data$expected_fits %>% 
    select(model,days_since_index,dow,expected=pred, expected_se=pred_se) %>% 
    left_join(select(data$excess_fits,-excess),by = join_by(model, days_since_index)) %>% 
    mutate_at(vars(excess_fit:excess_se), ~replace_na(.,0)) %>% 
    mutate(draw_mean = expected+excess_fit,
           draw_se = expected_se+excess_se)
  
  # ## Draw Error ----------------------------------------------------------------
  # if (empirical_error==TRUE){
  #   # Draw error term from empirical error values
  #   tmp_error <- map(data$prior_error_empirical, ~sample(.$error, size = data$upper_bound, replace = TRUE)) %>% 
  #     enframe(name = "model",value = "error") %>% 
  #     unnest("error") %>% 
  #     mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3)) 
  # } else {
  #   # Draw error term from distribution
  #   tmp_error <- data$prior_error_stats %>% 
  #     mutate(error = map(sd,~rnorm(n = data$upper_bound, mean = 0, sd =.))) %>% 
  #     unnest(error) %>% 
  #     select(model,error) %>% 
  #     mutate(days_since_index = rep((-data$upper_bound):(-1),times = 3))
  # }
  
  # # Draw Excess visits during delay window
  # tmp_excess <- data$excess_fits %>% 
  #   mutate(excess = map2_dbl(excess_fit,excess_se,~rnorm(1,mean = .x, sd = .y))) %>% 
  #   select(model,days_since_index, excess)
  # 
  # out <- tmp_expected %>% 
  #   inner_join(tmp_error, by = join_by(model, days_since_index)) %>% 
  #   left_join(tmp_excess, by = join_by(model, days_since_index)) %>% 
  #   mutate(excess = replace_na(excess, 0L)) %>% 
  #   mutate(simulated = expected+excess+error) %>% 
  #   select(model,days_since_index,dow,simulated)
  
  # Draw visits
  out <- tmp_expected %>% 
    mutate(simulated = map2_dbl(draw_mean,draw_se,~rnorm(1,mean = .x, sd = .y))) %>% 
    select(model,days_since_index,dow,simulated)
  
  return(out)
}

run_sim <- function(data, n_trials = 100){
  tibble(trial = 1:n_trials) %>% 
    mutate(simulated_data = map(trial,~simulate_visits(data = data)))
}


n_trials <- 1000

sim_res <- tibble(cp = cp_range) %>% 
  mutate(fits = map(cp,~find_fits(visit_counts,cp = .))) %>% 
  mutate(sim_data = map(fits,~run_sim(data = ., n_trials = n_trials)))



