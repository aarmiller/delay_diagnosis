library(tidyverse)
library(kableExtra)
library(ggplot2)
library(lubridate)


params <- list(proj = "measles")
# load delay_parms
load("/Volumes/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[params$proj]]
delay_params$cp <- 14+1
sim_in_path <- paste0("/Volumes/Statepi_Diagnosis/projects/", params$proj, "/delay_window_1_", delay_params$cp - 1,"/sim_results/")

load(paste0(sim_in_path,"fit_trends.RData"))

##################
#### Figure 1 ####
##################

fig1_A <- ssd_vis_count %>% 
  ggplot(aes(period,n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  scale_x_reverse() +
  theme_bw()+
  ylab("")+
  xlab("")

fig1_B <- ssd_vis_count %>% 
  ggplot(aes(period,n)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  geom_line(aes(y = pred), color = "red", linewidth = .7) +
  scale_x_reverse() +
  geom_vline(aes(xintercept = delay_params$cp-1), col = "blue", linewidth = .4) +
  theme_bw()+
  ylab("")+
  xlab("")

figure1 <- ggpubr::ggarrange(fig1_A, fig1_B,
                             labels = c("A", "B"),
                             ncol = 2, nrow = 1)
figure1 <- ggpubr::annotate_figure(figure1, left = grid::textGrob("Number of SSD Visits", rot = 90, vjust = 1, gp = grid::gpar(fontsize = 10)),
                bottom = grid::textGrob("Days Before Index Measles Diagnosis", gp = grid::gpar(fontsize = 10)),
                #top="Figure 1"
                ) 
ggsave("/Users/alanarakkal/Desktop/ongoing projects/measles/figure1.pdf",
       width = 5, height = 2.8,dpi = 1200,units = "in",
       plot = figure1)

# ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/measles/figures/figure1.pdf",
#        width = 5, height = 2.8,dpi = 600,units = "in",
#        plot = figure1)


# ggsave("publications/dengue/figures/figure1.jpg",
#        width = 5, height = 2.8,dpi = 600,units = "in",
#        plot = figure1)
# 
# ?ggsave


##################
#### Figure 2 ####
##################

ssd_vis_count %>% 
  filter(period<delay_params$cp) %>% 
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Measles Diagnosis") +
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10)) +
  # ggtitle("Figure 2",) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("/Users/alanarakkal/Desktop/ongoing projects/measles/figure2.pdf",
       width = 5, height = 2.8,dpi = 1200,units = "in")

# ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/measles/figures/figure2.pdf",
#        width = 5, height = 2.8, dpi = 600,units = "in")
# ?ggtitle
# 
# ggsave("publications/dengue/figures/figure2.jpg",
#        width = 5, height = 3,dpi = 600,units = "in")


########################
#### Temporal Trend ####
########################

time_data <- read_csv("/Volumes/Statepi_Diagnosis/projects/measles/delay_window_1_14/risk_models/month_year_cases.csv")
time_data2 <- read_csv("/Volumes/Statepi_Diagnosis/projects/measles/measles_validated/delay_window_1_14/risk_models/month_year_cases.csv")

time_data %>% 
  mutate(date = as.Date(paste(year, month, "1", sep = "-"))) %>% 
  ggplot(aes(x = date, y = n)) +
  geom_line(linewidth = 1) +
  geom_point( size = 2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    x = "Year",
    y = "Number of measles cases",
  ) +
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/alanarakkal/Desktop/ongoing projects/measles/figure3.pdf",
       width = 6, height = 2.8,dpi = 1200,units = "in")

cdc <- readxl::read_xlsx("/Users/alanarakkal/Desktop/ongoing projects/measles/data-table.xlsx")
time_data %>% 
  group_by(year) %>% 
  summarise(n = sum(n), .groups = "drop") %>% 
  mutate(Group = "Study Population") %>% 
  bind_rows(cdc %>% select(year, n = cases) %>% 
              mutate(Group = "CDC") ) %>% 
  bind_rows(time_data2 %>% 
              group_by(year) %>% 
              summarise(n = sum(n), .groups = "drop") %>% 
              mutate(Group = "Study Population - Validated") ) %>% 
  filter(between(year, 2001, 2023)) %>% 
  mutate(year = as.factor(year)) %>%
  ggplot(aes(x = year, y = n , colour = Group, group = Group)) +
  geom_line(linewidth = 1) +
  geom_point( size = 2) +
  labs(
    x = "Year",
    y = "Number of measles cases",
  ) +
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10)) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = NULL))

ggsave("/Users/alanarakkal/Desktop/ongoing projects/measles/figure3_v2.pdf",
       width = 6, height = 2.8,dpi = 1200,units = "in")

