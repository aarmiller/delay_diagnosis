library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)


cp <- 63+1

sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/blasto/delay_window_1_", cp - 1,"/sim_results/")
load(paste0(sim_in_path,"fit_trends.RData"))
ssd_vis_count_full <- ssd_vis_count


sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/blasto/blasto_inpatient/delay_window_1_", cp - 1,"/sim_results/")
load(paste0(sim_in_path,"fit_trends.RData"))
ssd_vis_count_inpatient <- ssd_vis_count


sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/blasto/blasto_outpatient/delay_window_1_", cp - 1,"/sim_results/")
load(paste0(sim_in_path,"fit_trends.RData"))
ssd_vis_count_outpatient <- ssd_vis_count


##################
#### Figure 1 ####
##################

fig1_A <- ssd_vis_count_full %>% mutate(Cohort = "Full Cohort") %>% 
  bind_rows(ssd_vis_count_inpatient %>% mutate(Cohort = "Inpatient Cohort"),
            ssd_vis_count_outpatient %>% mutate(Cohort = "Outpatient Cohort")) %>% 
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
  bind_rows(ssd_vis_count_inpatient %>% mutate(Cohort = "Inpatient Cohort"),
            ssd_vis_count_outpatient %>% mutate(Cohort = "Outpatient Cohort")) %>% 
  ggplot(aes(period,n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  geom_line(aes(y = pred), color = "red", linewidth = .7) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = cp-1), col = "blue", linewidth = .4) +
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
                bottom = grid::textGrob("Days Before Index Blastomycosis Diagnosis", gp = grid::gpar(fontsize = 11)),
                #top="Figure 1"
                ) 
ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/blasto/figures/figure1_combined.pdf",
       width = 6.5, height = 8,dpi = 600,units = "in",
       plot = figure1)

# ggsave("publications/blasto/figures/figure1.jpg",
#        width = 5, height = 2.8,dpi = 600,units = "in",
#        plot = figure1)
# 
# ?ggsave


##################
#### Figure 2 ####
##################

figure2 <- ssd_vis_count_full %>% mutate(Cohort = "Full Cohort") %>% 
  bind_rows(ssd_vis_count_inpatient %>% mutate(Cohort = "Inpatient Cohort"),
            ssd_vis_count_outpatient %>% mutate(Cohort = "Outpatient Cohort")) %>% 
  filter(period<cp) %>% 
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Blastomycosis Diagnosis") +
  theme_bw() +
  facet_wrap(vars(Cohort), ncol = 2, scales = "free_y") +
  scale_y_continuous(breaks = seq(0, 220, 20), limits = c(0, 220))+
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(face = "bold")) +
  # ggtitle("Figure 2",) +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/blasto/figures/figure2_combined.pdf",
       plot = figure2,
       width = 6, height = 5,dpi = 600,units = "in")

