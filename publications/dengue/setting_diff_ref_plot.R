library(tidyverse)

### Version 1 ------------------------------------------------------------------
load("/Volumes/Statepi_Diagnosis/projects/dengue/risk_models/reg_data/new_setting_labels_firth_g_2001_w_diff_ref_setting_summary_res.RData")
load("/Volumes/Statepi_Diagnosis/projects/dengue/risk_models/reg_data/new_setting_labels_firth_g_2001_w_ED_ref_setting.RData")

setting_labels <- expand.grid(outpatient = c("Outpatient", NA),
                              inpatient = c("Inpatient", NA),
                              ed = c("ED", NA),
                              obs_stay = c("Observational Stay", NA)) %>% 
  unite(., col = setting_label, sep = " + ", na.rm = T, remove =F) %>% 
  mutate(across(outpatient:obs_stay, ~ifelse(is.na(.), F, T))) %>% 
  mutate(setting_label = ifelse(setting_label == "", "Not any",
                                ifelse(setting_label == "Inpatient", "Inpatient Only",
                                       ifelse(setting_label == "Outpatient", "Outpatient Only",
                                              ifelse(setting_label == "ED", "ED Only",
                                                     ifelse(setting_label == "Observational Stay", "Observational Stay Only",
                                                            setting_label))))))
setting_labels <- setting_labels[order(rowSums(setting_labels[,2:5])), ]

temp <- c(setting_labels$setting_label[c(1:3, 5)], setting_labels %>% filter(grepl("Inpatient", setting_label)) %>% .$setting_label)
setting_factor <- c(temp, setting_labels$setting_label[!setting_labels$setting_label %in% temp])

plot_data <- outpatient_res %>% mutate(group = "Ref: Outpatient Only") %>% 
  bind_rows(inpatient_res %>% mutate(group = "Ref: Inpatient Only")) %>% 
  bind_rows(ed_res %>% mutate(group = "Ref: ED Only")) %>% 
  bind_rows(obs_res %>% mutate(group = "Ref: Observational Stay Only")) %>% 
  mutate(term = factor(term, levels = setting_factor),
         group = factor(group, levels = c("Ref: Observational Stay Only",
                                          "Ref: ED Only",
                                          "Ref: Outpatient Only",
                                          "Ref: Inpatient Only"))) 
# %>% 
#   mutate(across(est:high, ~ifelse(. == 1, NA, .)))


plot_data %>% 
  ggplot(aes(y = est, x = term))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust=1e-5, vjust = 0.5))+
  geom_hline(yintercept=1, linetype='dashed')+
  facet_grid(rows=vars(group), scales= "free_y")+
  ylab("Odds Ratio")+
  xlab("")
ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/dengue/figures/setting_diff_ref.pdf",
       height = 10,
       width = 8.2,
       units = "in")

### Version 2 ------------------------------------------------------------------
load("/Volumes/Statepi_Diagnosis/projects/dengue/risk_models/reg_data/new_setting_labels_firth_g_2001_w_diff_ref_setting_summary_res.RData")

outpatient_res <- outpatient_res %>% 
  mutate(term = stringr::str_replace_all(.$term, "Observational Stay", "Obs Stay"))

# Outpatient Plot
outpatient <- outpatient_res %>% 
  filter(grepl("Outpatient", term)) 

outpatient <- outpatient[c(7,1,6,2,3:5),] 

out_plot <- outpatient %>% 
  mutate(term = factor(term, levels = outpatient$term)) %>% 
  ggplot(aes(y = est, x = term))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust=1e-5, vjust = 0.5))+
  geom_hline(yintercept=1, linetype='dashed')+
  ylab("Odds Ratio")+
  xlab("")+
  ggtitle("Outpatient")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


temp_fun <- function(i, setting){
  temp <- trimws(stringr::str_split_1(i, "\\+"), "both")
  if(setting == "Inpatient" & temp[1] != setting){
    i <- paste0(c(temp[which(temp == setting)], temp[-which(temp == setting)]), collapse = " + ")
  } 
  if(setting != "Inpatient" & temp[1] != setting){
    i <- paste0(c(temp[which(temp == setting)],
                  temp[which(temp == "Inpatient")],
                  temp[-(which(temp == setting | temp == "Inpatient"))]), collapse = " + ")
  }
  return(i)
}

# Inpatient Plot

inpatient <- outpatient_res %>% 
  filter(grepl("Inpatient", term)) %>% 
  mutate(term2 = map_chr(term, ~temp_fun(i = .x, setting = "Inpatient"))) %>% 
  mutate(term2 = ifelse(term2 == "", term, term2))

inpatient <- inpatient[c(4,5,1,3,6,7,2),] %>% 
  select(term = term2, est, low, high)

in_plot <- inpatient %>% 
  mutate(term = factor(term, levels = inpatient$term)) %>% 
  ggplot(aes(y = est, x = term))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust=1e-5, vjust = 0.5))+
  geom_hline(yintercept=1, linetype='dashed')+
  ylab("Odds Ratio")+
  xlab("")+
  ggtitle("Inpatient")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# ED Plot

ED <- outpatient_res %>% 
  filter(grepl("ED", term)) %>% 
  mutate(term2 = map_chr(term, ~temp_fun(i = .x, setting = "ED"))) %>% 
  mutate(term2 = ifelse(term2 == "", term, term2))

ED <- ED[c(2,5,1,6,3,7,4),] %>% 
  select(term = term2, est, low, high)

ED_plot <- ED %>% 
  mutate(term = factor(term, levels = ED$term)) %>% 
  ggplot(aes(y = est, x = term))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust=1e-5, vjust = 0.5))+
  geom_hline(yintercept=1, linetype='dashed')+
  ylab("Odds Ratio")+
  xlab("")+
  ggtitle("ED")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Obs Plot

obs <- outpatient_res %>% 
  filter(grepl("Obs Stay", term)) %>% 
  mutate(term2 = map_chr(term, ~temp_fun(i = .x, setting = "Obs Stay"))) %>% 
  mutate(term2 = ifelse(term2 == "", term, term2))

obs <- obs[c(4,7,1,5,3,6,2),] %>% 
  select(term = term2, est, low, high)

obs_plot <- obs %>% 
  mutate(term = factor(term, levels = obs$term)) %>% 
  ggplot(aes(y = est, x = term))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust=1e-5, vjust = 0.5))+
  geom_hline(yintercept=1, linetype='dashed')+
  ylab("Odds Ratio")+
  xlab("")+
  ggtitle("Obs Stay")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

figure <- gridExtra::marrangeGrob(list(out_plot, ED_plot, obs_plot, in_plot),
                                  ncol = 2, nrow = 2,
                                  top = "")
ggsave("/Volumes/Statepi_Diagnosis/atlan/github/delay_diagnosis/publications/dengue/figures/setting_plot.png", 
       plot = figure, dpi = 500,
       height = 12, width = 12)




