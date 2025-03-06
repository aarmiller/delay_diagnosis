
library(tidyverse)
library(icd)
# Run Locally

rm(list = ls())

## Compute SSD Counts ----------------------------------------------------------

load("/Volumes/AML/tmp_transfer/jacob_pd_grant_all_dx_index.RData")

# Load SSD Codes
ssd_codes1 <- tibble(dx = c("3331", "7810", "7809", "7804", "7245", "78650", 
                           "29632", "78093", "7295", "32723", "311", "78039", 
                           "78009", "78057", "30000", "5990", "4019", "7802", 
                           "2989", "78052", "78791")) 

ssd_codes2 <- rename(read_csv("~/Data/jacob_pd_project/interesting_dx_ssds.csv"),dx=icd_9)

ssd_codes <- bind_rows(ssd_codes1,ssd_codes2) %>% 
  distinct() %>% 
  mutate(dx_ver = 9L) %>% 
  inner_join(codeBuildr::all_icd_labels) 


## Get Counts ------------------------------------------------------------------

obs_window <- 4*365

include_ids <- all_index %>% 
  filter(group == "PD") %>% 
  filter(time_before_index>=obs_window) %>% 
  distinct(patient_id)

ssd_counts <- all_dx_vis %>% 
  filter(group == "PD") %>% 
  filter(between(days_since_index,-obs_window, -1)) %>% 
  inner_join(include_ids, by = "patient_id") %>% 
  inner_join(ssd_codes, by = c("dx","dx_ver")) %>% 
  distinct(patient_id,dx,days_since_index) %>% 
  count(dx,days_since_index)

ssd_counts <- ssd_counts %>% 
  mutate(incidence = 1000*n/nrow(include_ids))

ssd_counts <- ssd_counts %>% 
  arrange(dx,days_since_index) %>% 
  group_by(dx) %>% 
  mutate(ma_n = (n + lag(n) + lag(n,2) +
                   lag(n,3) + lag(n,4) +
                   lag(n,5) + lag(n,6))/7) %>% 
  mutate(ma_incidence = 1000*n/nrow(include_ids)) %>% 
  ungroup()

ssd_counts_weekly <- ssd_counts %>% 
  select(dx,days_since_index,n) %>% 
  mutate(week = days_since_index %/% 7) %>% 
  group_by(dx,week) %>% 
  summarise(weekly_total = sum(n),
            weekly_incidence = 1000*weekly_total / nrow(include_ids)) %>% 
  ungroup()

ssd_counts_weekly <- ssd_counts_weekly %>% 
  filter(week> -209)

ssd_counts <- ssd_counts %>% 
  inner_join(ssd_codes) 

ssd_counts_weekly <- ssd_counts_weekly %>% 
  inner_join(ssd_codes)


#######################
#### Create Plots #####
#######################

#### Plot Export ---------------------------------------------------------------




plot_individual <- function(focus_ssd){
  focus_label <- ssd_codes %>% 
    filter(dx==focus_ssd) %>% 
    mutate(desc = paste0("ICD-9: ", dx,"- ",desc)) %>% 
    .$desc
  
  
  ssd_counts_weekly %>% 
    filter(dx == focus_ssd) %>% 
    ggplot(aes(-week,weekly_incidence)) +
    geom_line(size = 2) +
    theme_bw() +
    scale_x_reverse(breaks=c(52,104,156)) +
    xlab("Weeks Before Diagnosis") +
    ylab("Number of SSD Visits (per 1K patients)") +
    theme(axis.title = element_text(size = 14),
          axis.text  = element_text(size = 12)) +
    ggtitle(focus_label)
  ggsave(paste0("~/OneDrive - University of Iowa/grant_proposals/jacob_pd/figures/updated_ssd_plots/dx_",focus_ssd,".pdf"),
         width = 7, height = 5)
}

for(i in ssd_codes$dx){
  print(i)
  plot_individual(i)
}


#### Plot selected group -------------------------------------------------------

save(ssd_counts_weekly,file = "~/OneDrive - University of Iowa/grant_proposals/jacob_pd/data/ssd_counts_weekly.RData")



ssd_list <- c("3331","7810","7809","7804","7245","78650")

plot_data <- ssd_counts_weekly %>% 
  filter(dx %in% ssd_list) 

plot_data %>% 
  distinct(dx,desc)



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
