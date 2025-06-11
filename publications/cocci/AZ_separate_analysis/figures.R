
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)


load("/Volumes/AML/params/final_delay_params.RData")
all_data <- tibble()

for (i in c("cocci", "cocci_AZ", "cocci_not_AZ")){
# load delay_parms
params <- list(proj = i)

delay_params <- final_delay_params[[params$proj]]
delay_params$cp <- 91+1
sim_in_path <- paste0(delay_params$out_path, "delay_window_1_", delay_params$cp - 1,"/sim_results/")
sim_in_path <- str_replace(sim_in_path, "Shared", "Volumes")

load(paste0(sim_in_path,"fit_trends.RData"))

all_data <- bind_rows(all_data, 
                      ssd_vis_count %>% 
                        mutate(type = i)) 

rm(ssd_vis_count)

}

##################
#### Figure 1 ####
##################

figure1 <- all_data %>% 
  mutate(Population = factor(type, levels = c("cocci", "cocci_not_AZ", "cocci_AZ"),
                             labels = c("Full Cohort", "Non-AZ Cohort", "AZ Cohort"))) %>% 
  ggplot(aes(period, n , col = Population, group = Population)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp-1), col = "black", linewidth = .4) +
  theme_bw()+
  ylab("Number of SSD Visits")+
  xlab("Days Before Index Coccidiomycosis Diagnosis")

ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/cocci/AZ_separate_analysis/figures/figure1.pdf",
       width = 5, height = 2.8,dpi = 600,units = "in",
       plot = figure1)

figure2 <- all_data %>% 
  mutate(Population = factor(type, levels = c("cocci", "cocci_not_AZ", "cocci_AZ"),
                             labels = c("Full Cohort", "Non-AZ Cohort", "AZ Cohort"))) %>% 
  ggplot(aes(period, n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  geom_line(aes(y = pred), color = "red", linewidth = .7) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp-1), col = "black", linewidth = .4) +
  theme_bw()+
  facet_grid(row = vars(Population),
             scales = "free_y" )+
  ylab("Number of SSD Visits")+
  xlab("Days Before Index Coccidiomycosis Diagnosis")

ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/cocci/AZ_separate_analysis/figures/figure2.pdf",
       width = 5, height = 4,dpi = 600,units = "in",
       plot = figure2)
