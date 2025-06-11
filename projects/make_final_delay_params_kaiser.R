
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

final_delay_params <- list(sepsis_kaiser_cp_14= list(small_db_path = "/Shared/AML/kaiser_data/sepsis/delay_data/",
                                                     base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_kaiser/",  # base path to original prelim extract results
                                                     out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_kaiser/cp_14/",   # path to output delay new results
                                                     ssd_name = "sepsis_revised10",
                                                     cond_name = "sepsis_kaiser",
                                                     cp = 14+1,
                                                     upper_bound = 180,
                                                     final_model = "exponential",
                                                     periodicity = TRUE,
                                                     boot_trials = 100,
                                                     sim_trials = 100,
                                                     miss_bins = c(1,2,3,4,5),
                                                     duration_bins = 1:14),
                           
                           sepsis_kaiser_cp_7 = list(small_db_path = "/Shared/AML/kaiser_data/sepsis/delay_data/",
                                                     base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_kaiser/",  # base path to original prelim extract results
                                                     out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_kaiser/cp_7/",   # path to output delay new results
                                                     ssd_name = "sepsis_revised10",
                                                     cond_name = "sepsis_kaiser",
                                                     cp = 7+1,
                                                     upper_bound = 180,
                                                     final_model = "exponential",
                                                     periodicity = TRUE,
                                                     boot_trials = 100,
                                                     sim_trials = 100,
                                                     miss_bins = c(1,2,3,4,5),
                                                     duration_bins = 1:7)
                           
)

# save parameters
save(final_delay_params,file = "/Volumes/AML/params/final_delay_params_kaiser.RData")


