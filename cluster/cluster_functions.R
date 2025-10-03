plot_cluster_shapes <- function(cluster_data,k_val){
  tmp <- cluster_data$clust_res %>%
    filter(k == k_val)
  
  tmp$cluster_means[[1]] %>%
    ggplot(aes(days_since_index,mean_pred)) +
    geom_line(size = 1) +
    facet_wrap(~cluster) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    theme_minimal() +
    ylab("Normalized Visit Rate (Mean)") +
    xlab("Days Since Index") +
    ggtitle(paste0("Cluster Shapes for k = ",k_val))
}


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
4