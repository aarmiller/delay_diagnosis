
library(tidyverse)

make_abx_plots <- function(cond_name){
  
  ### Load delay params ----------------------------------------------------------
  
  load("/Shared/AML/params/delay_any_params.RData")
  
  delay_params <- delay_any_params[[cond_name]]
  
  base_path <- delay_params$path
  
  out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name)
  
  if (!dir.exists(out_path)){
    dir.create(out_path)
  }else{
    print("dir exists")
  }
  
  
  ### Paths for plot output ------------------------------------------------------
  
  plot_path <- paste0(out_path,"/abx_plots/")
  
  if (!dir.exists(plot_path)){
    dir.create(plot_path)
  }else{
    print("dir exists")
  }
  
  ### Load Data ------------------------------------------------------------------
  
  db <- src_sqlite(paste0(base_path,"/",cond_name,".db"))
  
  index_dx_dates <- db %>% tbl("index_dx_dates") %>% collect()
  all_rx_visits <- db %>% tbl("all_rx_visits") %>% collect()
  
  index_dx_dates <- index_dx_dates %>%
    filter(time_before_index>=delay_params$upper_bound)
  
  tmp_rx_visits <- all_rx_visits %>%
    inner_join(distinct(index_dx_dates, patient_id, index_date)) %>% 
    mutate(days_since_index = date - index_date) %>% 
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound))
  
  ## Get Antibiotics ---------------------------
  
  load("/Shared/AML/truven_mind_projects/antibiotic_risk_categories/antibiotics_groupings_new.RData")
  tmp_rx_visits <- tmp_rx_visits %>% 
    inner_join(antibiotic_ndc_groups_new)
  
  
  ### Make All abx plot -------------------------------------------------------
  tmp_rx_visits %>%
    distinct(patient_id,days_since_index) %>%
    count(days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) %>%
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    geom_smooth( color = "red",span = 0.25) +
    theme_minimal() +
    ggtitle(paste0("All abx before ",cond_name))
  ggsave(paste0(plot_path,"all_abx_before.pdf"),width = 6, height = 5)
  
  tmp_rx_visits %>%
    distinct(patient_id,days_since_index) %>%
    count(days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) %>%
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    geom_smooth( color = "red",span = 0.25) +
    theme_minimal() +
    ggtitle(paste0("All abx before and after ",cond_name))
  ggsave(paste0(plot_path,"all_abx_before_after.pdf"),width = 6, height = 5)
  
  ### Make abx by name plot -------------------------------------------------------
  by_name <- tmp_rx_visits %>%
    distinct(patient_id, name, days_since_index) %>%
    count(name, days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) 
  
  distinct_lab <- by_name %>% group_by(name) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% .$name
  split_lab <- split(distinct_lab, ceiling(seq_along(distinct_lab)/8))
  
  pdf(paste0(plot_path,"abx_by_name_before.pdf"),onefile = TRUE)
  
  for (i in 1:length(split_lab)){
    
    tmp_conds <- by_name %>%
      filter(name %in% split_lab[[i]])
    
    p <- tmp_conds %>%
      ggplot(aes(days_since_index,n), group = 1) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~name, scales = "free_y", nrow =4) +
      theme_minimal()
    
    print(p)
  }
  dev.off()
  
  by_name <- tmp_rx_visits %>%
    distinct(patient_id, name, days_since_index) %>%
    count(name, days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) 
  
  distinct_lab <- by_name %>% group_by(name) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% .$name
  split_lab <- split(distinct_lab, ceiling(seq_along(distinct_lab)/8))
  
  pdf(paste0(plot_path,"abx_by_name_before_after.pdf"),onefile = TRUE)
  
  for (i in 1:length(split_lab)){
    
    tmp_conds <- by_name %>%
      filter(name %in% split_lab[[i]])
    
    p <- tmp_conds %>%
      ggplot(aes(days_since_index,n), group = 1) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~name, scales = "free_y", nrow =4) +
      theme_minimal()
    
    print(p)
  }
  dev.off()
  
  
  ### Make abx by class plot -------------------------------------------------------
  by_class <- tmp_rx_visits %>%
    distinct(patient_id, class, days_since_index) %>%
    count(class, days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,-1)) 
  
  distinct_lab <- by_class %>% group_by(class) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% .$class
  split_lab <- split(distinct_lab, ceiling(seq_along(distinct_lab)/8))
  
  pdf(paste0(plot_path,"abx_by_class_before.pdf"),onefile = TRUE)
  
  for (i in 1:length(split_lab)){
    
    tmp_conds <- by_class %>%
      filter(class %in% split_lab[[i]])
    
    p <- tmp_conds %>%
      ggplot(aes(days_since_index,n), group = 1) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~class, scales = "free_y", nrow =4) +
      theme_minimal()
    
    print(p)
  }
  dev.off()
  
  by_class <- tmp_rx_visits %>%
    distinct(patient_id, class, days_since_index) %>%
    count(class, days_since_index) %>%
    filter(between(days_since_index,-delay_params$upper_bound,delay_params$upper_bound)) 
  
  distinct_lab <- by_class %>% group_by(class) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% .$class
  split_lab <- split(distinct_lab, ceiling(seq_along(distinct_lab)/8))
  
  pdf(paste0(plot_path,"abx_by_class_before_after.pdf"),onefile = TRUE)
  
  for (i in 1:length(split_lab)){
    
    tmp_conds <- by_class %>%
      filter(class %in% split_lab[[i]])
    
    p <- tmp_conds %>%
      ggplot(aes(days_since_index,n), group = 1) +
      geom_line() +
      geom_smooth( color = "red",span = 0.25) +
      facet_wrap(~class, scales = "free_y", nrow =4) +
      theme_minimal()
    
    print(p)
  }
  dev.off()
  
}