
library(tidyverse)
library(lubridate)
library(smallDB)
# devtools::install_github("aarmiller/smallDB")

print("started")

# name of condition
proj_name <- "cocci"
cond_name <- stringr::str_split(proj_name, "_")[[1]][1]

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[proj_name]]
delay_params$cp <- 91 + 1
delay_params$out_path <- paste0(final_delay_params[[proj_name]]$out_path, "delay_window_1_", delay_params$cp - 1, "/")

delay_base_path <- paste0(delay_params$base_path,"delay_results/")
sim_out_path <- paste0(delay_params$out_path,"sim_results/")

out_path <-paste0(delay_params$out_path,"risk_models/") 

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

## Load all cocci dx visits
load("/Shared/AML/truven_extracts/dx/cocci/cocci_dx_visits.RData")
load("/Shared/AML/truven_extracts/dx/cocci/cocci_enrolids.RData")

# connect db path
con <- DBI::dbConnect(RSQLite::SQLite(), 
                      paste0(delay_params$small_db_path, str_split(proj_name, "_")[[1]][1], ".db"))

# Get enrollment counts by state and year --------------------------------------
all_cocci_dx_visits <- all_outpatient_visits %>% distinct(enrolid, svcdate) %>% select(enrolid, admdate = svcdate) %>% 
  bind_rows(all_inpatient_visits %>% distinct(enrolid, admdate)) %>% 
  bind_rows(all_facility_visits %>% distinct(enrolid, svcdate) %>% select(enrolid, admdate = svcdate)) %>% 
  bind_rows(inpatient_services_visits %>% distinct(enrolid, admdate)) %>% 
  distinct(enrolid, admdate)

# 44,627 distinct enrollees all_cocci_dx_visits %>% distinct(enrolid), matches extract reports

# add year
all_cocci_dx_visits <- all_cocci_dx_visits %>% 
  mutate(year = as.character(year(as_date(admdate)))) 

# load crosswalk
load("/Shared/AML/truven_extracts/small_dbs/cocci/cocci_enrolid_crosswalk.RData")

all_cocci_dx_visits <- all_cocci_dx_visits %>% inner_join(enrolids) # only 44,292 unique enrolids in enrolids
# however there are supposed to be 44,627 enrollees diagnosed with cocci

# source(paste0(delay_params$out_path, "get_enroll_detail_fun.R"))
# load(paste0(delay_params$out_path, "egeoloc_labels.RData")) # checked with 2020 data dic on 07/31/2024
source(paste0(stringr::str_replace(delay_params$out_path,
                                   paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
              "get_enroll_detail_fun.R"))
load(paste0(stringr::str_replace(delay_params$out_path,
                                 paste0("delay_window_1_", delay_params$cp-1, "/"), ""),
            "egeoloc_labels.RData"))


enroll_collapsed_temp <- gather_collapse_enrollment(enrolid_list = all_cocci_dx_visits %>% distinct(patient_id) %>% .$patient_id,
                                                    vars = "egeoloc",
                                                    db_path =  paste0(delay_params$small_db_path,"cocci.db"),
                                                    num_cores=10,
                                                    collect_tab = collect_table(year = 1:22))

enroll_collapsed_temp2 <- enroll_collapsed_temp %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  inner_join(select(all_cocci_dx_visits,patient_id,admdate, year)) %>% 
  filter(admdate<=dtend & admdate>=dtstart) %>% 
  distinct() 
# only 43,209 out of 44,292 have location information, 1082 non-medicaid enrolids %>% count(medicaid) vs enroll_collapsed_temp2 %>% distinct(patient_id)
 
all_cocci_dx_visits <- all_cocci_dx_visits %>% 
  left_join(enroll_collapsed_temp2 %>% select(patient_id, egeoloc , admdate, year))

count_cocci_by_egeoloc_year <- all_cocci_dx_visits %>% 
  distinct() %>% 
  count(year, egeoloc) %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) 

## Now to compute truven enrollment by state and egeoloc -----------------------

get_egeoloc <- function(year){
  
  truven <- DBI::dbConnect(RSQLite::SQLite(), paste0("/Shared/Statepi_Marketscan/databases/Truven/truven_", year, ".db"))
  
  ccae <- tbl(truven, paste0("enrollment_detail_ccae_", year)) %>% select(enrolid, egeoloc) %>% 
    distinct(enrolid, egeoloc) %>% 
    count(egeoloc) %>% 
    collect(n = Inf)
  
  mdcr <- tbl(truven, paste0("enrollment_detail_mdcr_", year)) %>% select(enrolid, egeoloc) %>% 
    distinct(enrolid, egeoloc) %>% 
    count(egeoloc)%>% 
    collect(n = Inf)
  
  DBI::dbDisconnect(truven)

  out <- bind_rows(ccae, mdcr)  %>% 
    group_by(egeoloc) %>% 
    summarise(n = sum(n)) %>% 
    ungroup()
  
  return(out %>% mutate(year))
  
}

years_consider <- str_pad(1:22, width = 2, side = "left", pad = "0")

temp <- parallel::mclapply(years_consider,
                           mc.preschedule = F,
                           FUN = function(x){get_egeoloc(x)},
                           mc.cores = 22L)

truven_enrollement_by_egeoloc_year <- do.call("bind_rows", temp)

truven_cocci_inc <- truven_enrollement_by_egeoloc_year %>% 
  rename(n_total = n) %>% 
  mutate(year = paste0("20", year)) %>% 
  left_join(., count_cocci_by_egeoloc_year %>% rename(n_cases = n) %>% 
              select(year, n_cases, egeoloc), by = c("egeoloc", "year")) %>% 
  mutate(n_cases = ifelse(is.na(n_cases), 0, n_cases)) %>% 
  mutate(inc_per_100000 = (n_cases/n_total)*100000)

ave_truven_cocci_inc_state <- truven_cocci_inc %>% 
  filter(!is.na(egeoloc)) %>% 
  group_by(egeoloc) %>% 
  summarise(ave_inc_per_100000 = mean(inc_per_100000)) %>% 
  ungroup() %>% 
  inner_join(egeoloc_labels %>% select(egeoloc, location, state_name, state_abb)) %>% 
  arrange(desc(ave_inc_per_100000))

# save 
save(truven_cocci_inc, ave_truven_cocci_inc_state, file = paste0(out_path,"truven_inc_by_egeoloc_year.RData")) 

## Plot inc by state -----------------------------------------------------------
library(tidyverse)
library(usmap)
out_path <- "/Volumes/Statepi_Diagnosis/projects/cocci/delay_window_1_91/risk_models/"
load(paste0(out_path,"truven_inc_by_egeoloc_year.RData"))

ave_truven_cocci_inc_state <- ave_truven_cocci_inc_state %>% 
  mutate(state_abb=ifelse(location== "washington, dc" & is.na(state_abb), "DC", state_abb)) %>% 
  filter(!is.na(state_abb))

states_w_inc <- usmap::statepop %>% select(-pop_2015) %>%
  inner_join(ave_truven_cocci_inc_state %>% rename(abbr = state_abb))

plot_usmap(data = states_w_inc, values = "ave_inc_per_100000")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name =  "",
                       label = scales::comma) +
  labs(title = "Truven Coccidiomycosis Incidence per 100,000", subtitle = "Average annual incidence over 2001-2022") +
  theme(legend.position = "right")

# top 5
top5 <- states_w_inc %>% 
  arrange(desc(ave_inc_per_100000)) %>% 
  mutate(top5 = ifelse(row_number() <=5, 1,0))

plot_usmap(data = top5, values = "top5")+
  scale_fill_gradient2(low = "white",high = "green", name =  "",
                       label = scales::comma) +
  labs(title = "Top 5 - Truven Coccidiomycosis Incidence per 100,000", subtitle = "Average annual incidence over 2001-2022") +
  theme(legend.position = "none")

# top 10
top10 <- states_w_inc %>% 
  arrange(desc(ave_inc_per_100000)) %>% 
  mutate(top10 = ifelse(row_number() <=10, 1,0))

plot_usmap(data = top10, values = "top10")+
  scale_fill_gradient2(low = "white",high = "green", name =  "",
                       label = scales::comma) +
  labs(title = "Top 10 - Truven Coccidiomycosis Incidence per 100,000", subtitle = "Average annual incidence over 2001-2022") +
  theme(legend.position = "none")

# boxplot(states_w_inc$ave_inc_per_100000)

