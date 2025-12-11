rm(list = ls())
library(tidyverse)

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/cluster_results/2025_12_01/cluster_res.RData")

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/cluster_results/dx_counts.RData")

dx_counts

visit_counts


isolate_cluster <- function(cluster_data,k_val,cluster_number,dx_ver){
  tmp <- cluster_data$clust_res %>%
    filter(k == k_val) %>%
    select(clust_labels) %>%
    unnest(clust_labels) %>%
    filter(cluster==cluster_number) %>%
    rename(rank = index) %>%
    inner_join(distinct(cluster_data$pred_data,rank,code),
               by = "rank")
  
  if (dx_ver == 9){
    tmp %>%
      inner_join(filter(codeBuildr::all_icd_labels,dx_ver==9) %>%
                   select(code = dx, desc))
  } else {
    tmp %>%
      inner_join(filter(codeBuildr::all_icd_labels,dx_ver==10) %>%
                   select(code = dx, desc))
  }
  
}


library(icd)

find_focal_clusters <- function(cluster_data,focal_code,dx_ver){
  
  focal_index <- distinct(cluster_data$pred_data,code,index) %>%
    filter(code == focal_code)
  
  focal_clusters <- cluster_data$clust_res %>%
    select(k,clust_labels) %>%
    unnest(clust_labels) %>%
    filter(index == focal_index$index) %>%
    distinct(k,cluster)
  
  tmp_means <- cluster_data$clust_res %>%
    select(k,cluster_means) %>%
    unnest(cluster_means) %>%
    inner_join(focal_clusters,by = join_by(k, cluster))
  
  tmp_conds <-  focal_clusters %>%
    inner_join(cluster_data$clust_res %>%
                 select(k,clust_labels) %>%
                 unnest(clust_labels),
               by = join_by(k, cluster)) %>%
    inner_join(distinct(cluster_data$pred_data,code,index),
               by = join_by(index)) %>%
    count(code,rank=index) %>%
    arrange(desc(n),rank)
  
  if (dx_ver == 9){
    tmp_conds <- tmp_conds %>%
      inner_join(filter(codeBuildr::all_icd_labels,dx_ver==9) %>%
                   select(code = dx, desc),
                 by = join_by(code))
  } else {
    tmp_conds <- tmp_conds %>%
      inner_join(filter(codeBuildr::all_icd_labels,dx_ver==10) %>%
                   select(code = dx, desc),
                 by = join_by(code))
  }
  
  list(cluster_means = tmp_means,
       condition_counts = tmp_conds)
}

#### Cough ####

cough_10 <- find_focal_clusters(cluster_data = clust_res_10,
                               focal_code = "R05",
                               dx_ver = 10L)


cough_9 <- find_focal_clusters(cluster_data = clust_res_9,
                               focal_code = "7862",
                               dx_ver = 9L)

#### Shortness of breath ####

sob_10 <- find_focal_clusters(cluster_data = clust_res_10,
                              focal_code = "R0602",
                              dx_ver = 10L)

sob_9 <- find_focal_clusters(cluster_data = clust_res_9,
                              focal_code = "78605",
                              dx_ver = 9L)


#### Weight Loss ####

weight_loss_10 <- find_focal_clusters(cluster_data = clust_res_10,
                                      focal_code = "R634",
                                      dx_ver = 10L)

weight_loss_9 <- find_focal_clusters(cluster_data = clust_res_9,
                                      focal_code = "78321",
                                      dx_ver = 9L)






focal_cluster_res_10 <- bind_rows(mutate(cough_10$condition_counts,term = "cough"),
                                  mutate(sob_10$condition_counts,term = "shortness_of_breath"),
                                  mutate(weight_loss_10$condition_counts,term = "weight_loss"))


focal_cluster_res_10 %>% 
  group_by(code,rank,desc) %>% 
  summarise(total = sum(n)) %>% 
  arrange(desc(total),rank) %>% 
  ungroup() %>% 
  mutate(cluster_frac = round(100*total/(3*28),2)) %>% 
  select(code,desc,`Overall Rank`=rank,`Total Clusters`=total,`Cluster Fraction`=cluster_frac) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_res_icd_10.csv")

focal_cluster_res_9 <- bind_rows(mutate(cough_9$condition_counts,term = "cough"),
                                  mutate(sob_9$condition_counts,term = "shortness_of_breath"),
                                  mutate(weight_loss_9$condition_counts,term = "weight_loss"))

focal_cluster_res_9 %>% 
  group_by(code,rank,desc) %>% 
  summarise(total = sum(n)) %>% 
  arrange(desc(total),rank) %>% 
  ungroup() %>% 
  mutate(cluster_frac = round(100*total/(3*28),2)) %>% 
  select(code,desc,`Overall Rank`=rank,`Total Clusters`=total,`Cluster Fraction`=cluster_frac) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_res_icd_9.csv")




### Plot cluster Conditions ###

top_conds9 <- focal_cluster_res_9 %>% 
  group_by(code,rank,desc) %>% 
  summarise(total = sum(n)) %>% 
  arrange(desc(total),rank) %>% 
  ungroup() %>% 
  mutate(cluster_frac = round(100*total/(3*28),2))

top_conds10 <- focal_cluster_res_10 %>% 
  group_by(code,rank,desc) %>% 
  summarise(total = sum(n)) %>% 
  arrange(desc(total),rank) %>% 
  ungroup() %>% 
  mutate(cluster_frac = round(100*total/(3*28),2))

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/cluster_results/2025_12_01/cluster_count_data.RData")


tmp_dx9 <- dx9_counts_nested %>% 
  select(code,data) %>% 
  unnest(data) %>% 
  inner_join(icd9_ranks) %>% 
  ungroup() %>% 
  mutate(label = paste0(code,": ",str_sub(description,1,40)))

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_plots_icd9.pdf",onefile = TRUE)
for (i in 1:ceiling(nrow(top_conds9)/6)){
  tmp_conds <- top_conds9 %>%
    slice((i*6-5):(i*6)) %>% 
    select(code)
  
  p <- tmp_conds %>%
    inner_join(tmp_dx9, by ="code") %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,frac, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~label, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none") +
    ylab("Percentage of Visits") +
    xlab("Days Since Index")
  
  print(p)
}
dev.off()


### Plot dx 10  before and after -----------------------------------------------


tmp_dx10 <- dx10_counts_nested %>% 
  select(code,data) %>% 
  unnest(data) %>% 
  inner_join(icd10_ranks) %>% 
  ungroup() %>% 
  mutate(label = paste0(code,": ",str_sub(description,1,40)))

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_plots_icd10.pdf",onefile = TRUE)
for (i in 1:ceiling(2000/6)){
  tmp_conds <- icd10_ranks %>%
    slice((i*6-5):(i*6)) %>% 
    select(code)
  
  p <- tmp_conds %>%
    inner_join(tmp_dx10, by ="code") %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,frac, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~label, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none") +
    ylab("Percentage of Visits") +
    xlab("Days Since Index")
  
  print(p)
}
dev.off()

#### Export excluded conditions ####

excluded_codes9 <- dx9_counts_nested %>% 
  select(code) %>% 
  ungroup() %>% 
  anti_join(top_conds9) 

excluded_codes9 %>% 
  inner_join(filter(codeBuildr::all_icd_labels, dx_ver==9) %>% 
               select(code = dx, desc)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_res_excluded_icd9.csv")

top_conds9

excluded_codes10 <- dx10_counts_nested %>% 
  select(code) %>% 
  ungroup() %>% 
  anti_join(top_conds10) 

excluded_codes10 %>% 
  inner_join(filter(codeBuildr::all_icd_labels, dx_ver==10) %>% 
               select(code = dx, desc)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_res_excluded_icd10.csv")


pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_excluded_plots_icd9.pdf",onefile = TRUE)
for (i in 1:ceiling(nrow(excluded_codes9)/6)){
  tmp_conds <- excluded_codes9 %>%
    slice((i*6-5):(i*6)) %>% 
    select(code)
  
  p <- tmp_conds %>%
    inner_join(tmp_dx9, by ="code") %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,frac, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~label, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none") +
    ylab("Percentage of Visits") +
    xlab("Days Since Index")
  
  print(p)
}
dev.off()

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/cluster_excluded_plots_icd10.pdf",onefile = TRUE)
for (i in 1:ceiling(nrow(excluded_codes10)/6)){
  tmp_conds <- excluded_codes10 %>%
    slice((i*6-5):(i*6)) %>% 
    select(code)
  
  p <- tmp_conds %>%
    inner_join(tmp_dx10, by ="code") %>%
    filter(days_since_index!=0) %>%
    mutate(period = ifelse(days_since_index<0,"Before","After")) %>%
    ggplot(aes(days_since_index,frac, color = period)) +
    geom_line() +
    geom_smooth(aes(group = period),color = "red",span = 0.25) +
    facet_wrap(~label, scales = "free_y", nrow =4) +
    geom_vline(aes(xintercept = 0)) +
    theme_minimal() +
    theme(legend.position="none") +
    ylab("Percentage of Visits") +
    xlab("Days Since Index")
  
  print(p)
}
dev.off()




