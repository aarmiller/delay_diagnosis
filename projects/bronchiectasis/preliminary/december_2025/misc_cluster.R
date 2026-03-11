rm(list = ls())
library(tidyverse)

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/cluster_results/2025_12_29/1095_1095/cluster_res.RData")

load("~/Data/Statepi_Diagnosis/prelim_results/bronchiectasis/cluster_results/2025_12_29/1095_1095/cluster_count_data.RData")

base_path <- paste0("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/",
                    str_replace_all(Sys.Date(),"-","_"))

if(!dir.exists(base_path)){dir.create(base_path)}

out_path <- paste0(base_path,"/set7/")

if(!dir.exists(out_path)){dir.create(out_path)}


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

cough_10 <- find_focal_clusters(cluster_data = clust_res_10_pre,
                                focal_code = "R05",
                                dx_ver = 10L)


cough_9 <- find_focal_clusters(cluster_data = clust_res_9_pre,
                               focal_code = "7862",
                               dx_ver = 9L)

#### Shortness of breath ####

sob_10 <- find_focal_clusters(cluster_data = clust_res_10_pre,
                              focal_code = "R0602",
                              dx_ver = 10L)

sob_9 <- find_focal_clusters(cluster_data = clust_res_9_pre,
                             focal_code = "78605",
                             dx_ver = 9L)


#### Chest Pain ####

chest_pain_10 <- find_focal_clusters(cluster_data = clust_res_10_pre,
                                     focal_code = "R079",
                                     dx_ver = 10L)


chest_pain_9 <- find_focal_clusters(cluster_data = clust_res_9_pre,
                                    focal_code = "78650",
                                    dx_ver = 9L)

clust_res_9_pre$clust_res$clust_labels[[1]]
filter(code == "4940")

#### Hemoptysis ####

hemoptysis_10 <- find_focal_clusters(cluster_data = clust_res_10_pre,
                                     focal_code = "R042",
                                     dx_ver = 10L)


hemoptysis_9 <- find_focal_clusters(cluster_data = clust_res_9_pre,
                                    focal_code = "7863",
                                    dx_ver = 9L)

focal_counts9 <- bind_rows(mutate(hemoptysis_9$condition_counts,focal = "hemoptysis"),
          mutate(chest_pain_9$condition_counts,focal = "chest_pain"),
          mutate(sob_9$condition_counts,focal = "sob"),
          mutate(cough_9$condition_counts,focal = "cough"))


focal_counts10 <- bind_rows(mutate(hemoptysis_10$condition_counts,focal = "hemoptysis"),
                           mutate(chest_pain_10$condition_counts,focal = "chest_pain"),
                           mutate(sob_10$condition_counts,focal = "sob"),
                           mutate(cough_10$condition_counts,focal = "cough"))

focal_counts9 %>% 
  group_by(code,desc) %>% 
  summarise(rank = mean(rank),
            max_n = max(n/28),
            mean_n = mean(n/28),
            min_n = min(n/28)) %>% 
  arrange(desc(max_n),desc(mean_n),rank) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/2026_01_06/rank_include9.csv")

focal_counts10 %>% 
  group_by(code,desc) %>% 
  summarise(rank = mean(rank),
            max_n = max(n/28),
            mean_n = mean(n/28),
            min_n = min(n/28)) %>% 
  arrange(desc(max_n),desc(mean_n),rank) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/bronchiectasis/prelim_results/2026_01_06/rank_include10.csv")
