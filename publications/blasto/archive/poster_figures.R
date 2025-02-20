

library(tidyverse)

cond_name <- "blasto"

sim_out_path <- "/Volumes/Statepi_Diagnosis/prelim_results/blasto/delay_results/ssd_visit/"


#### SSD Trend plot ####


load("/Volumes/Statepi_Diagnosis/prelim_results/blasto/delay_results/ssd_visit/ssd_fit_res.RData")

ssd_fit_res$model_fits %>% 
  filter(model=="quadratic") %>% 
  ggplot(aes(period,n)) +
  geom_line() +
  theme_bw() +
  scale_x_reverse() +
  geom_vline(aes(xintercept = 60), linetype = 2) +
  geom_line(aes(y = value), color = "red") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))


#### Histogram ####

ssd_fit_res$model_fits %>% 
  # filter(model=="quadratic") %>% 
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))

