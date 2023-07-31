
################################################################################
### This Script Contains the parameters necessary for running automated builds 
### and analysis 
################################################################################


##################################
#### Preliminary Delay Params ####
##################################

### Example --------------------------------------------------------------------
# list(path = "/Shared/AML/small_dbs/<disease>/truven/enroll_restrict_365/",
#      cp = 30,
#      cp_lower = NA,                  # OPTIONAL: Specify the lowest change-point to evaluate
#      cp_upper = NA,                  # OPTIONAL: Specify the highest change-point to evaluate
#      upper_bound = 365,
#      final_model = NA,               # options include c("linear", "quadratic", "cubic", "exponential")
#      periodicity = TRUE,
#      miss_bins = c(1,2,3,4,5),
#      duration_bins = c(1,2,3,4,5,6,7,10,14,17))


delay_any_params <- list(cvst = list(path = "/Shared/AML/truven_extracts/small_dbs/cvst/",
                                      cp = 21,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 180,
                                      final_model = NA,
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                         
                         histo = list(path = "/Shared/AML/truven_extracts/small_dbs/histo/",
                                      cp = 100,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 365,
                                      final_model = NA,
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60)),
                         
                         meningitis = list(path = "/Shared/AML/truven_extracts/small_dbs/meningitis/",
                                           cp = 50,
                                           cp_lower = NA,
                                           cp_upper = NA,
                                           upper_bound = 365,
                                           final_model = NA,
                                           periodicity = TRUE,
                                           miss_bins = c(1,2,3,4,5),
                                           duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60)),
                         
                         sepsis_revised10 = list(path = "/Shared/AML/truven_extracts/small_dbs/sepsis_revised10/",
                                                 cp = 14,
                                                 cp_lower = NA,
                                                 cp_upper = NA,
                                                 upper_bound = 180,
                                                 final_model = NA,
                                                 periodicity = TRUE,
                                                 miss_bins = c(1,2,3,4,5),
                                                 duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                         
                         sarcoid = list(path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                      cp = 180,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 365*2,
                                      final_model = NA,
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90,120,150,180)))

# save parameters
save(delay_any_params,file = "/Volumes/AML/params/delay_any_params.RData")

# update master list
library(tidyverse)
load("/Volumes/AML/params/delay_any_params.RData")

length(delay_any_params)
tmp <- names(delay_any_params)

# create folders for holding output
for (i in 1:length(delay_any_params)){
  if(!dir.exists(paste0("/Volumes/Statepi_Diagnosis/prelim_results/",tmp[i]))){
    dir.create(paste0("/Volumes/Statepi_Diagnosis/prelim_results/",tmp[i]))
  }
}

