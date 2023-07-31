
rm(list = ls())
library(tidyverse)


load(paste0("/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/delay_results/all_dx_visits.RData"))

cond_name <- "meningitis_bacterial"

load("/Shared/AML/params/final_delay_params.RData")

delay_params <- final_delay_params[[cond_name]]

rm(final_delay_params)


### Specific 


#### OLD (from sepsis remove)

# filter to specific index cases
# 1) 180 days of continuous enrollment
# 2) all prior data in the ICD-10 era (180 days after Jan 1 2016)

# index_cases <- index_dx_dates %>% 
#   filter(time_before_index>=delay_params$upper_bound) %>% 
#   filter(index_date>=as.integer(ymd("2016-01-01")))
# 
# save(index_cases, file = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/index_cases.RData")

