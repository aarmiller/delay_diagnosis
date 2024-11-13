

rm(list = ls())
library(tidyverse)

cond_name <- "sepsis_pre_covid"

load("/Volumes/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

sim_out_path <- paste0("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/","sim_results/")

# Load Index cases
load("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/index_cases.RData")

n_patients <- nrow(index_cases)


tmp <- read_csv("/Volumes/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/sim_results/exponential_cp14/model1_estimate_num_before.csv")

tmp %>% 
  mutate_at(vars(or:conf.high),~round(.,2)) %>% 
  mutate(ci = paste0(conf.low,"-",conf.high)) %>% 
  select(term, or, ci) %>% 
  mutate(term = str_remove(term,"TRUE")) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/results/reg_results.csv") 


#### Trend plot ####


load(paste0(sim_out_path,"report_data/trends_output.RData"))
load(paste0(sim_out_path,"trends/fit_trends.RData"))

trends_ssd <- ssd_vis_count_fitted %>% 
  rename(tot_miss=num_miss) %>% 
  unnest(counts) %>% 
  mutate(model_label = paste0(model," (CP = ",cp," days)"))

# mse_res_ssd <- trends_ssd %>%
#   filter(is.na(num_miss)) %>% 
#   group_by(model_label) %>% 
#   summarise(rmse = sqrt(mean((pred-n)^2))) %>% 
#   mutate(label = paste("RMSE: ", round(rmse,2)))
# 
# trends_all <- all_vis_count_fitted %>% 
#   rename(tot_miss=num_miss) %>% 
#   unnest(counts) %>% 
#   mutate(model_label = paste0(model," (CP = ",cp," days)"))
# 
# mse_res_all <- trends_all %>%
#   filter(is.na(num_miss)) %>% 
#   group_by(model_label) %>% 
#   summarise(rmse = sqrt(mean((pred-n)^2))) %>% 
#   mutate(label = paste("RMSE: ", round(rmse,2)))
# 

### Appendix Figure 1 ----------------------------------------------------------

y_pos <- .9*max(trends_ssd$n)
x_pos <- .8*180

mse_res_ssd <- trends_ssd %>%
  filter(is.na(num_miss)) %>%
  mutate(model_label = ifelse(cp==7,"7-day opportunity window",
                              "14-day opportunity window")) %>% 
  group_by(model_label) %>%
  summarise(rmse = sqrt(mean((pred-n)^2))) %>%
  mutate(label = paste("RMSE: ", round(rmse,2)))

mse_res_ssd <- mse_res_ssd %>% 
  mutate(key = NA)


trends_ssd %>%
  select(period,"Observed"=n,"Expected"=pred,cp) %>% 
  gather(key = key,value = value,-period,-cp) %>% 
  mutate(model_label = ifelse(cp==7,"7-day opportunity window",
                              "14-day opportunity window")) %>% 
  ggplot(aes(period,value)) +
  geom_line(aes(color = key)) +
  scale_x_reverse() +
  facet_wrap(~model_label) +
  geom_vline(aes(xintercept = cp), linetype =2) +
  theme_bw() +
  scale_color_manual(values = c("red","black")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  geom_text(data = mse_res_ssd,
            mapping = aes(x = x_pos, y = y_pos, label = label)) +
  ylab("Number of SSD visits") +
  xlab("Days before index sepsis diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/figures/appendix_figure1.pdf",
       width = 8, height = 5)

 ### Figure 1 -------------------------------------------------------------------

trends_ssd %>%
  filter(cp ==14) %>%  
  select(period,"Observed"=n,"Expected"=pred,cp) %>% 
  gather(key = key,value = value,-period,-cp) %>% 
  ggplot(aes(period,value,color = key)) +
  geom_line() +
  scale_x_reverse() +
  geom_vline(aes(xintercept = cp), linetype =2) +
  theme_bw() +
  scale_color_manual(values = c("red","black")) +
  ylab("Total visits with signs or symptoms") +
  xlab("Days before index sepsis diagnosis") +
  theme(legend.title = element_blank(),
        legend.position = c(0.25, 0.8))
# ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/figures/figure1.pdf",
#        width = 4.5, height = 4)
ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/submissions/crit_care_explorations/figure1.pdf",
       width = 4.5, height = 4, dpi = 600)
# 
# ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/submissions/crit_care_explorations/figure1.tiff",dpi = 300,
#        width = 4.5, height = 4)


### Methods Graphic (Revisions) ------------------------------------------------

trends_ssd %>%
  filter(cp ==14) %>%  
  filter(period<=75) %>% 
  select(period,"Observed"=n,"Expected"=pred,cp)  %>% 
  ggplot(aes(period,Observed)) +
  geom_point(size = 2) +
  geom_vline(aes(xintercept = cp), linetype =2) +
  # geom_vline(aes(xintercept = 1), linetype =1) +
  theme_bw() +
  geom_ribbon(aes(ymin=0,ymax = Expected-425), fill = "blue", alpha = 0.5) +
  geom_ribbon(aes(ymin=Expected,ymax = Observed), fill = "red", alpha = 0.5) +
  scale_color_manual(values = c("red","black")) +
  geom_line(aes(y = Expected), color = "blue", size = 1.2) +
  ylab("Total SSD visits") +
  xlab("Days Before Index Sepsis Diagnosis") +
  scale_x_reverse(breaks = c(1,14,75),expand = c(0.014, 0)) +
  scale_y_continuous(expand = c(0.021,0))
ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/submissions/crit_care_explorations/revisions/methods_ssd_fig.pdf",
       width = 4.7, height = 3, dpi = 600)

ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/submissions/crit_care_explorations/revisions/methods_ssd_fig.jpg",
       width = 4.7, height = 3, dpi = 600)


### Figure 2 -------------------------------------------------------------------



trends_ssd %>% 
  filter(cp ==14) %>%  
  ggplot(aes(period,num_miss)) +
  geom_histogram(stat = "identity") +
  scale_x_reverse() +
  ylab("Number of Missed Opportunities") +
  xlab("Days Before Index Diagnosis") +
  theme_bw() 
ggsave("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/submissions/crit_care_explorations/figure2.pdf",
       width = 4.5, height = 4, dpi = 600)
