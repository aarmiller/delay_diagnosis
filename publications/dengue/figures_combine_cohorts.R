library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)


params <- list(proj = "dengue")
# load delay_parms
load("/Volumes/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[params$proj]]
delay_params$cp <- 14+1
sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/dengue/delay_window_1_", delay_params$cp - 1,"/sim_results/")

load(paste0(sim_in_path,"fit_trends.RData"))

ssd_vis_count_full <- ssd_vis_count


params <- list(proj = "dengue_validated")
# load delay_parms
load("/Volumes/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[params$proj]]
delay_params$cp <- 14+1
sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/dengue/dengue_validated/delay_window_1_", delay_params$cp - 1,"/sim_results/")

load(paste0(sim_in_path,"fit_trends.RData"))

ssd_vis_count_validated <- ssd_vis_count

##################
#### Figure 1 ####
##################

fig1_A <- ssd_vis_count_full %>% mutate(Cohort = "Full Cohort") %>% 
  bind_rows(ssd_vis_count_validated %>% mutate(Cohort = "Validated Cohort")) %>% 
  ggplot(aes(period,n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  scale_x_reverse() +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free_y") +
  theme_bw()+
  theme(strip.text= element_text(face = "bold"),
        axis.title.y = element_text(size = 10))+
  ylab("Number of SSD Visits")+
  xlab("")


fig1_B <- ssd_vis_count_full %>% mutate(Cohort = "Full Cohort") %>% 
  bind_rows(ssd_vis_count_validated %>% mutate(Cohort = "Validated Cohort")) %>% 
  ggplot(aes(period,n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  geom_line(aes(y = pred), color = "red", linewidth = .7) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp-1), col = "blue", linewidth = .4) +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free_y") +
  theme_bw()+
  theme(strip.text = element_text(face = "bold"),
        axis.title.y = element_text(size = 10))+
  ylab("Number of SSD Visits")+
  xlab("")


figure1 <- ggpubr::ggarrange(fig1_A, fig1_B,
                             label.x = 0,
                             label.y = 1,
                             labels = c("A", "B"),
                             ncol = 1, nrow = 2)
figure1 <- ggpubr::annotate_figure(figure1,,
                bottom = grid::textGrob("Days Before Index Dengue Diagnosis", gp = grid::gpar(fontsize = 11)),
                #top="Figure 1"
                ) 
ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/dengue/figures/figure1_combined.pdf",
       width = 6, height = 5,dpi = 600,units = "in",
       plot = figure1)

# ggsave("publications/dengue/figures/figure1.jpg",
#        width = 5, height = 2.8,dpi = 600,units = "in",
#        plot = figure1)
# 
# ?ggsave


##################
#### Figure 2 ####
##################

figure2 <- ssd_vis_count_full %>% mutate(Cohort = "Full Cohort") %>% 
  bind_rows(ssd_vis_count_validated %>% mutate(Cohort = "Validated Cohort")) %>% 
  filter(period<delay_params$cp) %>% 
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Dengue Diagnosis") +
  theme_bw() +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free_y") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(face = "bold")) +
  # ggtitle("Figure 2",) +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/dengue/figures/figure2_combined.pdf",
       plot = figure2,
       width = 6, height = 5,dpi = 600,units = "in")

