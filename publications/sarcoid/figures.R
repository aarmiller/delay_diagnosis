library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)


params <- list(proj = "sarcoid")
cond_name <- "Sarcoid"
# load delay_parms
load("/Volumes/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[params$proj]]
sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/", params$proj, "/sim_results/")


load(paste0(sim_in_path,"fit_trends.RData"))

##################
#### Figure 1 ####
##################

fig1_A <- ssd_vis_count %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line() +
  scale_x_reverse() +
  theme_bw()+
  ylab("")+
  xlab("")

fig1_B <- ssd_vis_count %>% 
  ggplot(aes(period,n)) +
  geom_point() +
  geom_line() +
  geom_line(aes(y = pred), color = "red", linewidth = 1.1) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp), col = "blue", linewidth = 1.0) +
  theme_bw()+
  ylab("")+
  xlab("")

figure1 <- ggpubr::ggarrange(fig1_A, fig1_B,
                             labels = c("A", "B"),
                             ncol = 2, nrow = 1)
figure1 <- ggpubr::annotate_figure(figure1, left = grid::textGrob("Number of SSD Visits", rot = 90, vjust = 1, gp = grid::gpar(fontsize = 10)),
                bottom = grid::textGrob(paste0("Days Before Index ", cond_name, " Diagnosis"), gp = grid::gpar(fontsize = 10))) 
ggsave(paste0("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", params$proj, "/figures/figure1.pdf"),
       width = 10, height = 6,
       plot = figure1)

##################
#### Figure 2 ####
##################

ssd_vis_count %>% 
  filter(period<delay_params$cp) %>% 
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab(paste0("Days Before Index ", cond_name, " Diagnosis")) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12))
ggsave(paste0("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/", params$proj, "/figures/figure2.pdf"),
       width = 10, height = 6)

