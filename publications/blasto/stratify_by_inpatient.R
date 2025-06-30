

db <- src_sqlite("~/Data/MarketScan/truven_extracts/small_dbs/blasto/blasto.db")

tmp <- codeBuildr::load_disease_codes("blasto",return_tibble = TRUE)

blasto_dx <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(tmp, copy = TRUE) %>% 
  collect()

blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date))


blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date)) %>% 
  inner_join(blasto_dx) %>% 
  filter(index_dx==date) %>% 
  filter(inpatient==1)


blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date)) %>% 
  inner_join(blasto_dx) %>% 
  filter(date>=index_dx & date<=(index_dx+30)) %>% 
  filter(inpatient==1)


ssd_codes <- codeBuildr::load_ssd_codes("blasto") %>% 
  mutate(dx_ver = ifelse(type == "icd9", 9L, 10L)) %>% 
  select(dx = code, dx_ver)


ssd_vis <- db %>% 
  tbl("all_dx_visits") %>% 
  inner_join(ssd_codes, copy = TRUE) %>% 
  filter(between(days_since_index,-365,-1)) %>% 
  collect()



include_ids <- db %>% 
  tbl("index_dx_dates") %>% 
  collect() %>% 
  filter(time_before_index>=365)

outpatient_ssd_days <- ssd_vis %>% 
  filter(inpatient==0) %>% 
  distinct(patient_id,days_since_index)

admission_window <- 30

admitted_cases <- blasto_dx %>% 
  group_by(patient_id) %>% 
  summarise(index_dx = min(date)) %>% 
  inner_join(blasto_dx) %>% 
  filter(date>=index_dx & date<=(index_dx+admission_window)) %>% 
  filter(inpatient==1) %>% 
  distinct(patient_id) %>% 
  mutate(inpatient=1L)


cases <- select(include_ids,patient_id,index_date) %>% 
  left_join(admitted_cases) %>% 
  mutate(inpatient = replace_na(inpatient,0L)) %>% 
  select(patient_id,inpatient)

outpatient_ssd_days %>% 
  inner_join(cases) %>% 
  group_by(inpatient) %>% 
  count(days_since_index) %>% 
  ungroup() %>% 
  inner_join(count(cases,inpatient, name = "tot_pop")) %>% 
  mutate(frac = n/tot_pop) %>% 
  mutate(`Case Type` = ifelse(inpatient == 1,"Admitted Case","Outpatient Case")) %>% 
  ggplot(aes(-days_since_index,frac, color = `Case Type`)) +
  geom_line() +
  scale_x_reverse() +
  theme_bw() +
  xlab("Days Before Blasto Diagnosis") +
  ylab("Percent of Cases with outpatient SSD Visit") +
  ggtitle("Outpatient SSD visits before blasto diagnosis",
          subtitle = paste0(admission_window,"-day inpatient admission window from diagnosis"))


count_cases <- function(admission_window){
  admitted_cases <- blasto_dx %>% 
    group_by(patient_id) %>% 
    summarise(index_dx = min(date)) %>% 
    inner_join(blasto_dx) %>% 
    filter(date>=index_dx & date<=(index_dx+admission_window)) %>% 
    filter(inpatient==1) %>% 
    distinct(patient_id) %>% 
    mutate(inpatient=1L)
  
  
  cases <- select(include_ids,patient_id,index_date) %>% 
    left_join(admitted_cases) %>% 
    mutate(inpatient = replace_na(inpatient,0L)) %>% 
    select(patient_id,inpatient)
  
  outpatient_ssd_days %>% 
    inner_join(cases) %>% 
    group_by(inpatient) %>% 
    filter(days_since_index<=63) %>% 
    summarise(n_patients = n_distinct(patient_id),
              n_visits = n()) %>% 
    inner_join(count(cases,inpatient, name = "tot_pop")) %>% 
    mutate(frac_patients = round(100*n_patients/tot_pop,2),
           mean_visits = round(n_visits/tot_pop,2)) %>% 
    mutate(out1 = paste0(n_patients," (",frac_patients,")"),
           out2 = paste0(n_visits, " (",mean_visits,")")) %>% 
    select(inpatient,tot_pop,out1,out2)
  
}



case_counts <- tibble(window = c(0,7,14,21,30)) %>% 
  mutate(res = map(window,count_cases)) %>% 
  unnest(res)


case_counts %>% 
  write_csv("~/Downloads/case_counts.csv")
