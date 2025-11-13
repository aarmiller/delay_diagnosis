rm(list = ls())
library(tidyverse)

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/cocci/cocci.db")

index_dates <- db %>% 
  tbl("index_dx_dates") %>% 
  collect()


# ids with sufficient enrollment
include_ids <- index_dates %>% 
  filter(time_before_index>=365)

abx_codes <- codeBuildr::load_rx_codes("all_abx") %>% 
  enframe() %>% 
  unnest(value) %>% 
  rename(ndc = value)


collect_codes <- distinct(abx_codes, ndcnum = ndc)

abx_vis <- db %>% 
  tbl("all_rx_visits") %>% 
  inner_join(collect_codes, copy = TRUE) %>% 
  collect()


abx_counts <- abx_vis %>% 
  inner_join(rename(abx_codes,ndcnum=ndc)) %>% 
  distinct(patient_id,date,name) %>% 
  inner_join(select(include_ids,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,name,days_since_index) %>% 
  count(name,days_since_index) 

abx_counts <- distinct(abx_counts,name) %>% 
  mutate(days_since_index = map(name,~-365:-1)) %>% 
  unnest(days_since_index) %>% 
  left_join(abx_counts) %>% 
  mutate(n = replace_na(n,0L))

abx_totals <- abx_counts %>% 
  group_by(name) %>% 
  summarise(abx_total = sum(n)) %>% 
  arrange(desc(abx_total))

#### Any ABX --------------

abx_vis %>% 
  inner_join(select(include_ids,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index) %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  theme_bw() +
  xlab("Days before Cocci Diagnosis") +
  ylab("Number of Patients prescribed Antibiotic")


abx_vis %>% 
  inner_join(select(include_ids,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  mutate(week = (days_since_index %/% 7)) %>% 
  filter(week> -53) %>% 
  count(week) %>% 
  ggplot(aes(week,n)) +
  geom_line() +
  theme_bw() +
  xlab("Weeks before Cocci Diagnosis") +
  ylab("Number of Patients prescribed Antibiotic")

### Final Counts ###

abx_trends_daily <- abx_vis %>% 
  inner_join(select(include_ids,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  count(days_since_index)

abx_trends_weekly <- abx_vis %>% 
  inner_join(select(include_ids,patient_id,index_date)) %>% 
  mutate(days_since_index = date-index_date) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  mutate(week = (days_since_index %/% 7)) %>% 
  filter(week> -53) %>% 
  count(week)


### Compare to SSD Visits ###

ssd_codes <- codeBuildr::load_ssd_codes("cocci") %>% 
  mutate(dx_ver = ifelse(type == "icd9",9L,10L)) %>% 
  select(dx = code,dx_ver)


ssd_counts <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(ssd_codes,copy = TRUE, by = join_by(dx, dx_ver)) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  distinct(patient_id,days_since_index) %>% 
  collect() %>% 
  count(days_since_index,name = "ssd_visits")

ssd_counts %>% 
  ggplot(aes(days_since_index,ssd_visits)) +
  geom_line()

bind_rows(mutate(abx_trends_daily,Group = "Antibiotic Prescriptions"),
          rename(ssd_counts,n = ssd_visits) %>% 
            mutate(Group = "SSD Visits")) %>% 
  group_by(Group) %>% 
  mutate(norm_n = (n-min(n))/(max(n)-min(n))) %>% 
  ggplot(aes(-days_since_index,norm_n, color = Group)) +
  geom_line(size = 1) +
  theme_bw() +
  scale_x_reverse() +
  xlab("Days Before Coccidioidomycosis Diagnosis") +
  ylab("Number of Visits Normalized") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
ggsave("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/submissions/OFID/figure3.pdf",
       width = 12, height = 10,units = "in",dpi = 300)

abx_trends_daily %>% 
  ggplot(aes(-days_since_index,n)) +
  geom_line(size = 1) +
  theme_bw() +
  scale_x_reverse() +
  xlab("Days Before Coccidioidomycosis Diagnosis") +
  ylab("Number of Patients with Antibiotic Prescription") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
ggsave("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/submissions/OFID/appendix_fig2.pdf",
         width = 12, height = 8,units = "in", dpi = 300)



#### By Type --------------


nrow(abx_totals)/6


i <- 1


abx_counts %>% 
  inner_join(tmp_conds, by = "name") %>% 
  ggplot(aes(days_since_index,n)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", nrow = 3) +
  theme_bw()

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/results/abx_before.pdf",onefile = TRUE)
for (i in 1:12){
  tmp_conds <- abx_totals %>%
    slice((i*6-5):(i*6)) %>% 
    select(name)
  
  p <- tmp_conds %>% 
    inner_join(abx_counts) %>% 
    ggplot(aes(days_since_index,n)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    theme_bw()
  
  print(p)
}
dev.off()

pdf("~/OneDrive - University of Iowa/WorkingPapers/delay_dx_projects/cocci/results/abx_before_weekly.pdf",onefile = TRUE)
for (i in 1:12){
  tmp_conds <- abx_totals %>%
    slice((i*6-5):(i*6)) %>% 
    select(name)
  
  p <- tmp_conds %>% 
    inner_join(abx_counts,by = join_by(name)) %>% 
    mutate(week = -(days_since_index %/% 7)) %>% 
    filter(week<53) %>% 
    group_by(name,week) %>% 
    summarise(n = sum(n)) %>% 
    ungroup() %>% 
    ggplot(aes(-week,n)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y", nrow = 3) +
    theme_bw()
  
  print(p)
}
dev.off()
