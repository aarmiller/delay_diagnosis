
# return fitted values and number missed for a given function
return_fits <- function(data, model, cp, periodicity = TRUE){
  
  if (periodicity==TRUE){
    
    if (model == "linear"){
      
      fit <- lm(n ~ period + dow, filter(data, period>cp))
      
    } else if (model == "exponential") {
      
      fit <- lm(n ~ log(period) + dow, filter(data, period>cp))
      
    } else if (model == "quadratic") {
      
      fit <- lm(n ~ poly(period,degree = 2) + dow, filter(data, period>cp))
      
    } else if (model == "cubic") {
      
      fit <- lm(n ~ poly(period,degree = 3) + dow, filter(data, period>cp))
      
    }
    
  } else {
    
    if (model == "linear"){
      
      fit <- lm(n ~ period, filter(data, period>cp))
      
    } else if (model == "exponential") {
      
      fit <- lm(n ~ log(period), filter(data, period>cp))
      
    } else if (model == "quadratic") {
      
      fit <- lm(n ~ poly(period,degree = 2), filter(data, period>cp))
      
    } else if (model == "cubic") {
      
      fit <- lm(n ~ poly(period,degree = 3), filter(data, period>cp))
      
    }
    
  }
  
  data %>%
    mutate(pred = predict(fit,newdata = .)) %>% 
    mutate(num_miss = n-pred) %>%
    mutate(num_miss = ifelse(num_miss<0,0,num_miss)) %>%
    mutate(num_miss = ifelse(period>cp,NA,num_miss))
  
}