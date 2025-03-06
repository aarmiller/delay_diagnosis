
library(tidyverse)
load("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/data/ssd_counts_weekly.RData")
     
# list of ssds to plot
ssd_list <- c("3331","7810","7809","7804","7245","78650")

# plot data for selected ssds
plot_data <- ssd_counts_weekly %>% 
  filter(dx %in% ssd_list) 

# plot
plot_data %>% 
  ggplot(aes(-week,weekly_incidence)) +
  geom_line(size = 1) +
  facet_wrap(~desc, scale = "free_y") +
  theme_bw() +
  scale_x_reverse(breaks=c(52,104,156)) +
  xlab("Weeks Before Diagnosis") +
  ylab("Number of SSD Visits (per 1K patients)") +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        strip.text = element_text(size = 12))
ggsave(paste0("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/figures/selected_ssds.pdf"),
       width = 12, height = 7)
ggsave(paste0("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/figures/selected_ssds.svg"),
       width = 12, height = 7)
