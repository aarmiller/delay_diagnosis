

rm(list = ls())
library(tidyverse)
library(bit64)
library(parallel)
library(smallDB)
library(delaySim)
library(delayCluster)

# devtools::install_github("aarmiller/delayCluster", dependencies = FALSE, force = TRUE)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "ami"

# set number of dx codes to consider
dx_cluter_lim <- 1000

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

out_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name)

cluster_out_path <- paste0(out_path,"/cluster_results")


# load cluster data
load(paste0(cluster_out_path,"/dx_counts.RData"))

###########################
#### Correct DX Counts ####
###########################

dx_count_bounds <- summarise(dx_counts,max_bound = max(days_since_index),
                             min_bound =min(days_since_index))

## Correct ICD-9 Counts --------------------------------------------------------

dx9_counts <- dx_counts %>%
  filter(dx_ver==9) %>%
  rename(code=dx) %>%
  select(-dx_ver)

# Create count holders for missing values
dx9_span <- dx9_counts %>%
  distinct(code) %>%
  mutate(tmp = map(code, ~tibble(days_since_index=dx_count_bounds$min_bound:dx_count_bounds$max_bound))) %>%
  unnest(tmp)

dx9_counts <- dx9_span %>%
  left_join(dx9_counts) %>%
  select(-n_visits) %>%
  left_join(filter(visit_counts,dx_ver==9)) %>%
  select(code,days_since_index,n_dx,n_visits=n,frac) %>%
  mutate(n_dx = replace_na(n_dx,0),
         frac = replace_na(frac,0))

## Correct ICD-10 Counts -------------------------------------------------------

dx10_counts <- dx_counts %>%
  filter(dx_ver==10) %>%
  rename(code=dx) %>%
  select(-dx_ver)

# Create count holders for missing values
dx10_span <- dx10_counts %>%
  distinct(code) %>%
  mutate(tmp = map(code, ~tibble(days_since_index=dx_count_bounds$min_bound:dx_count_bounds$max_bound))) %>%
  unnest(tmp)

dx10_counts <- dx10_span %>%
  left_join(dx10_counts) %>%
  select(-n_visits) %>%
  left_join(filter(visit_counts,dx_ver==10)) %>%
  select(code,days_since_index,n_dx,n_visits=n,frac) %>%
  mutate(n_dx = replace_na(n_dx,0),
         frac = replace_na(frac,0))


#############################
#### Get Condition Ranks ####
#############################

tmp <- dx_counts %>%
  filter(between(days_since_index,-14,-1)) %>%
  group_by(dx,dx_ver) %>%
  summarise(n = sum(n_dx)) %>%
  arrange(dx_ver,desc(n)) %>%
  ungroup()

icd9_ranks <- tmp %>%
  filter(dx_ver==9) %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number()) %>%
  left_join(icd::icd9cm_hierarchy %>% as_tibble() %>%
              select(dx=code,long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  select(code = dx, n_14day =n, rank, description =long_desc)

icd10_ranks <- tmp %>%
  filter(dx_ver==10) %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number()) %>%
  left_join(codeBuildr::icd10cm_labels %>% as_tibble() %>%
              select(dx=code,long_desc) %>%
              mutate(dx = as.character(dx))) %>%
  select(code = dx, n_14day =n, rank, description =long_desc)



##########################
#### Fit Loess Models ####
##########################


### fit loess models icd-9 -----------------------------------------------------

dx9_counts_nested <- icd9_ranks %>%
  filter(rank<=dx_cluter_lim) %>%
  select(code) %>%
  inner_join(dx9_counts) %>%
  filter(days_since_index!=0) %>%
  group_by(code) %>%
  mutate(norm_n = (frac-min(frac))/(max(frac)-min(frac))) %>%
  ungroup() %>%
  select(code,frac,norm_n,days_since_index) %>%
  group_by(code) %>%
  nest()

dx9_counts_nested <- dx9_counts_nested %>%
  mutate(data1 = map(data,~filter(.,days_since_index<0))) %>%
  mutate(data2 = map(data,~filter(.,days_since_index>0))) %>%
  mutate(fit_pre = map(data1,~loess(norm_n~days_since_index, data = .,span = 0.2))) %>%
  mutate(fit_post = map(data2,~loess(norm_n~days_since_index, data = .,span = 0.2))) %>%
  mutate(data1 = map2(data1,fit_pre,~mutate(.x,pred = predict(.y)))) %>%
  mutate(data2 = map2(data2,fit_post,~mutate(.x,pred = predict(.y))))

dx9_fits <- bind_rows(dx9_counts_nested %>%
                        select(code,data1) %>%
                        unnest(data1) %>%
                        ungroup(),
                      dx9_counts_nested %>%
                        select(code,data2) %>%
                        unnest(data2) %>%
                        ungroup()) %>%
  arrange(code,days_since_index) %>%
  inner_join(select(icd9_ranks,code,index = rank,description))

### fit loess models icd-10 -----------------------------------------------------

dx10_counts_nested <- icd10_ranks %>%
  filter(rank<=1000) %>%
  select(code) %>%
  inner_join(dx10_counts) %>%
  filter(days_since_index!=0) %>%
  group_by(code) %>%
  mutate(norm_n = (frac-min(frac))/(max(frac)-min(frac))) %>%
  ungroup() %>%
  select(code,frac,norm_n,days_since_index) %>%
  group_by(code) %>%
  nest()

dx10_counts_nested <- dx10_counts_nested %>%
  mutate(data1 = map(data,~filter(.,days_since_index<0))) %>%
  mutate(data2 = map(data,~filter(.,days_since_index>0))) %>%
  mutate(fit_pre = map(data1,~loess(norm_n~days_since_index, data = .,span = 0.2))) %>%
  mutate(fit_post = map(data2,~loess(norm_n~days_since_index, data = .,span = 0.2))) %>%
  mutate(data1 = map2(data1,fit_pre,~mutate(.x,pred = predict(.y)))) %>%
  mutate(data2 = map2(data2,fit_post,~mutate(.x,pred = predict(.y))))

dx10_fits <- bind_rows(dx10_counts_nested %>%
                         select(code,data1) %>%
                         unnest(data1) %>%
                         ungroup(),
                       dx10_counts_nested %>%
                         select(code,data2) %>%
                         unnest(data2) %>%
                         ungroup()) %>%
  arrange(code,days_since_index) %>%
  inner_join(select(icd10_ranks,code,index = rank,description))

########################
#### Find Clusters #####
########################

#### functions -----------------------------------------------------------------
find_pred_clusters <- function(pred_data, k_range=3:30){

  loess_fits_trans <- pred_data %>%
    select(index,code,description,days_since_index,pred) %>%
    mutate(days_since_index = paste0("X",days_since_index)) %>%
    spread(key = days_since_index,value = pred)

  clust_res <- tibble(k=k_range) %>%
    mutate(clusters = map(k,~find_clusters(fit_data = loess_fits_trans, k =.x)))

  clust_res <- clust_res %>%
    mutate(clust_labels = map(clusters,~.$distances %>% select(index,cluster) %>% ungroup()),
           wss = map_dbl(clusters,~.$wss))

  clust_means <- clust_res %>%
    select(k, clust_labels) %>%
    unnest(clust_labels) %>%
    inner_join(pred_data, by = "index") %>%
    group_by(k,cluster,days_since_index) %>%
    summarise(mean_pred = mean(pred)) %>%
    group_by(k) %>%
    nest() %>%
    rename(cluster_means = data)

  clust_res <- clust_res %>%
    inner_join(clust_means, by = "k") %>%
    select(k,clust_labels:cluster_means)

  out <- list(clust_res = clust_res,
              pred_data = pred_data)

  return(out)


}

### ICD-9 Clusters -------------------------------------------------------------

clust_res_9 <- find_pred_clusters(dx9_fits, k_range = 3:30)


### ICD-9 Clusters -------------------------------------------------------------

clust_res_10 <- find_pred_clusters(dx10_fits, k_range = 3:30)


save(clust_res_9,clust_res_10,file = paste0(cluster_out_path,"/cluster_res.RData"))
