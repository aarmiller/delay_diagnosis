
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


delay_any_params <- list(ami = list(path = "/Shared/AML/truven_extracts/small_dbs/ami/",
                                       cp = 30,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 365,
                                       final_model = "cubic",
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
                         blasto = list(path = "/Shared/AML/truven_extracts/small_dbs/blasto/",
                                     cp = 51,
                                     cp_lower = NA,
                                     cp_upper = NA,
                                     upper_bound = 365,
                                     final_model = "quadratic",
                                     periodicity = TRUE,
                                     miss_bins = c(1,2,3,4,5),
                                     duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
                         chf = list(path = "/Shared/AML/truven_extracts/small_dbs/chf/",
                                       cp = 90,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 365*2,
                                       final_model = "quadratic",
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
                         cocci = list(path = "/Shared/AML/truven_extracts/small_dbs/cocci/",
                                       cp = 60,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 365,
                                       final_model = "cubic",
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60)),
                         
                         cvst = list(path = "/Shared/AML/truven_extracts/small_dbs/cvst/",
                                      cp = 21,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 180,
                                      final_model = NA,
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                         
                         dengue = list(path = "/Shared/AML/truven_extracts/small_dbs/dengue/",
                                       cp = 30,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 180,
                                       final_model = NA,
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                         
                         endocarditis = list(path = "/Shared/AML/truven_extracts/small_dbs/endocarditis/",
                                       cp = 30,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 180,
                                       final_model = NA,
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                         
                         histo = list(path = "/Shared/AML/truven_extracts/small_dbs/histo/",
                                      cp = 70,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 365,
                                      final_model = "cubic",
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60)),
                         
                         lung_cancer = list(path = "/Shared/AML/truven_extracts/small_dbs/lung_cancer/",
                                   cp = 120,
                                   cp_lower = NA,
                                   cp_upper = NA,
                                   upper_bound = 2*365,
                                   final_model = "cubic",
                                   periodicity = TRUE,
                                   miss_bins = c(1,2,3,4,5),
                                   duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
                         meningitis = list(path = "/Shared/AML/truven_extracts/small_dbs/meningitis/",
                                           cp = 50,
                                           cp_lower = NA,
                                           cp_upper = NA,
                                           upper_bound = 365,
                                           final_model = NA,
                                           periodicity = TRUE,
                                           miss_bins = c(1,2,3,4,5),
                                           duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60)),
                         
                         pcp = list(path = "/Shared/AML/truven_extracts/small_dbs/pcp/",
                                       cp = 75,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 365,
                                       final_model = NA,
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
                         pe = list(path = "/Shared/AML/truven_extracts/small_dbs/pe/",
                                       cp = 98,
                                       cp_lower = NA,
                                       cp_upper = NA,
                                       upper_bound = 365,
                                       final_model = "quadratic",
                                       periodicity = TRUE,
                                       miss_bins = c(1,2,3,4,5),
                                       duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,50)),
                         
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
                                      cp = 128,
                                      cp_lower = NA,
                                      cp_upper = NA,
                                      upper_bound = 365*2,
                                      final_model = "quadratic",
                                      periodicity = TRUE,
                                      miss_bins = c(1,2,3,4,5),
                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90,120,150,180)),
                         
                         tb = list(path = "/Shared/AML/truven_extracts/small_dbs/tb/",
                                        cp = 90,
                                        cp_lower = NA,
                                        cp_upper = NA,
                                        upper_bound = 365,
                                        final_model = "quadratic",
                                        periodicity = TRUE,
                                        miss_bins = c(1,2,3,4,5),
                                        duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90,120,150,180)),
                         
                         pertussis = list(path = "/Shared/AML/truven_extracts/small_dbs/pertussis/",
                                        cp = 21,
                                        cp_lower = NA,
                                        cp_upper = NA,
                                        upper_bound = 365,
                                        final_model = NA,
                                        periodicity = TRUE,
                                        miss_bins = c(1,2,3,4,5),
                                        duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90,120,150,180))
                         )

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

