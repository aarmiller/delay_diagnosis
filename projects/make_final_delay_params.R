
################################################################################
### This Script Contains the parameters necessary for the final project analysis 
################################################################################

# Notes: Set this up to work with two types of sub analysis:
#        1) Subset of patients (e.g., Resp TB vs other TB)
#        2) Shifted index dates (e.g., evidence of disease before diagnosis)
#
# For the above to work we need to create subsets of the datasets used in the 
# preliminary simulations. The script to generate these recreated object should 
# be stored in delay_diagnosis/projects/<disease>
  
final_delay_params <- list(sarcoid = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                          base_path = "/Shared/Statepi_Diagnosis/prelim_results/sarcoid/",  # base path to original prelim extract results
                                          out_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/",   # path to output delay new results
                                          cp = 150,
                                          upper_bound = 365*2,
                                          final_model = NA,
                                          periodicity = TRUE,
                                          boot_trials = 100,
                                          sim_trials = 100,
                                          miss_bins = c(1,2,3,4,5,7,10),
                                          duration_bins = c(3,7,14,21,30,60,90,120,150)),
                           
                           sarcoid_lung = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                               base_path = "/Shared/Statepi_Diagnosis/prelim_results/sarcoid/",  # base path to original prelim extract results
                                               out_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/",   # path to output delay new results
                                               cp = 180,
                                               upper_bound = 365*2,
                                               final_model = "cubic",
                                               periodicity = TRUE,
                                               boot_trials = 100,
                                               sim_trials = 100,
                                               miss_bins = c(1,2,3,4,5,7,10),
                                               duration_bins = c(3,7,14,21,30,60,90,120,150)),
                           
                           sarcoid_skin = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                               base_path = "/Shared/Statepi_Diagnosis/prelim_results/sarcoid/",  # base path to original prelim extract results
                                               out_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_skin/",   # path to output delay new results
                                               cp = 180,
                                               upper_bound = 365*2,
                                               final_model = "cubic",
                                               periodicity = TRUE,
                                               boot_trials = 100,
                                               sim_trials = 100,
                                               miss_bins = c(1,2,3,4,5,7,10),
                                               duration_bins = c(3,7,14,21,30,60,90,120,150))
                           
)

# save parameters
save(final_delay_params,file = "/Volumes/AML/params/final_delay_params.RData")


