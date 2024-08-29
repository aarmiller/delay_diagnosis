
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
  
final_delay_params <- list(meningitis_bacterial = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/meningitis/",
                                                       base_path = "/Shared/Statepi_Diagnosis/prelim_results/meningitis/",  # base path to original prelim extract results
                                                       out_path = "/Shared/Statepi_Diagnosis/projects/meningitis_bacterial/",   # path to output delay new results
                                                       ssd_name = "meningitis",
                                                       cp = 50,
                                                       upper_bound = 365,
                                                       final_model = NA,
                                                       periodicity = TRUE,
                                                       boot_trials = 100,
                                                       sim_trials = 100,
                                                       miss_bins = c(1,2,3,4,5,7,10),
                                                       duration_bins = c(3,7,14,21,30,40,50)),
                           
                           sarcoid = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                          base_path = "/Shared/Statepi_Diagnosis/prelim_results/sarcoid/",  # base path to original prelim extract results
                                          out_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/",   # path to output delay new results
                                          ssd_name = "sarcoid",
                                          cp = 128,
                                          upper_bound = 365*2,
                                          final_model = "cubic",
                                          periodicity = TRUE,
                                          boot_trials = 100,
                                          sim_trials = 100,
                                          miss_bins = c(1,2,3,4,5),
                                          duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90,120,150,180)),
                           
                           sarcoid_lung = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sarcoid/",
                                               base_path = "/Shared/Statepi_Diagnosis/prelim_results/sarcoid/",  # base path to original prelim extract results
                                               out_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/",   # path to output delay new results
                                               ssd_name = "sarcoid_lung",
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
                                               ssd_name = "sarcoid_skin",
                                               cp = 180,
                                               upper_bound = 365*2,
                                               final_model = "cubic",
                                               periodicity = TRUE,
                                               boot_trials = 100,
                                               sim_trials = 100,
                                               miss_bins = c(1,2,3,4,5,7,10),
                                               duration_bins = c(3,7,14,21,30,60,90,120,150)),
                           
                           sepsis_revised10 = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sepsis_revised10/",
                                                   base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/",  # base path to original prelim extract results
                                                   out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/",   # path to output delay new results
                                                   ssd_name = "sepsis_revised10",
                                                   cp = c(7,14),
                                                   upper_bound = 180,
                                                   final_model = c("exponential","cubic"),
                                                   periodicity = TRUE,
                                                   boot_trials = 100,
                                                   sim_trials = 100,
                                                   miss_bins = c(1,2,3,4,5),
                                                   duration_bins = c(1:14)),
                           
                           sepsis_pre_covid = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sepsis_revised10/",
                                                   base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/",  # base path to original prelim extract results
                                                   out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/pre_covid/",   # path to output delay new results
                                                   ssd_name = "sepsis_revised10",
                                                   cp = c(7,14),
                                                   upper_bound = 180,
                                                   final_model = "exponential",
                                                   periodicity = TRUE,
                                                   boot_trials = 100,
                                                   sim_trials = 100,
                                                   miss_bins = c(1,2,3,4,5),
                                                   duration_bins = c(1:14)),
                           
                           sepsis_early_covid = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sepsis_revised10/",
                                                     base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/",  # base path to original prelim extract results
                                                     out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/early_covid/",   # path to output delay new results
                                                     ssd_name = "sepsis_revised10",
                                                     cp = c(7,14),
                                                     upper_bound = 180,
                                                     final_model = "exponential",
                                                     periodicity = TRUE,
                                                     boot_trials = 100,
                                                     sim_trials = 100,
                                                     miss_bins = c(1,2,3,4,5),
                                                     duration_bins = c(1:14)),
                           
                           sepsis_covid = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/sepsis_revised10/",
                                               base_path = "/Shared/Statepi_Diagnosis/prelim_results/sepsis_revised10/",  # base path to original prelim extract results
                                               out_path = "/Shared/Statepi_Diagnosis/projects/sepsis_revised10/covid/",   # path to output delay new results
                                               ssd_name = "sepsis_revised10",
                                               cp = c(7,14),
                                               upper_bound = 180,
                                               final_model = "exponential",
                                               periodicity = TRUE,
                                               boot_trials = 100,
                                               sim_trials = 100,
                                               miss_bins = c(1,2,3,4,5),
                                               duration_bins = c(1:14)),
                           
                           dengue = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/dengue/",
                                         base_path = "/Shared/Statepi_Diagnosis/prelim_results/dengue/",  # base path to original prelim extract results
                                         out_path = "/Shared/Statepi_Diagnosis/projects/dengue/",   # path to output delay new results
                                         ssd_name = "dengue",
                                         cp = 17,
                                         upper_bound = 180,
                                         final_model = "cubic",
                                         periodicity = TRUE,
                                         boot_trials = 100,
                                         sim_trials = 100,
                                         miss_bins = c(1,2,3,4,5),
                                         duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                           
                           
                           dengue_validated = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/dengue/",
                                                   base_path = "/Shared/Statepi_Diagnosis/prelim_results/dengue/",  # base path to original prelim extract results
                                                   out_path = "/Shared/Statepi_Diagnosis/projects/dengue/dengue_validated/",   # path to output delay new results
                                                   ssd_name = "dengue",
                                                   cp = 17,
                                                   upper_bound = 180,
                                                   final_model = "cubic",
                                                   periodicity = TRUE,
                                                   boot_trials = 100,
                                                   sim_trials = 100,
                                                   miss_bins = c(1,2,3,4,5),
                                                   duration_bins = c(1,2,3,4,5,6,7,10,14,17,21)),
                           
                           
                           blasto = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/blasto/",
                                         base_path = "/Shared/Statepi_Diagnosis/prelim_results/blasto/",  # base path to original prelim extract results
                                         out_path = "/Shared/Statepi_Diagnosis/projects/blasto/",   # path to output delay new results
                                         ssd_name = "blasto",
                                         cp = 92,
                                         upper_bound = 365,
                                         final_model = "cubic",
                                         periodicity = TRUE,
                                         boot_trials = 100,
                                         sim_trials = 100,
                                         miss_bins = c(1,2,3,4,5),
                                         duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90)),
                           
                           blasto_top2_baddley = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/blasto/",
                                                      base_path = "/Shared/Statepi_Diagnosis/prelim_results/blasto/",  # base path to original prelim extract results
                                                      out_path = "/Shared/Statepi_Diagnosis/projects/blasto/blasto_top2_baddley/",   # path to output delay new results
                                                      ssd_name = "blasto",
                                                      cp = 92,
                                                      upper_bound = 365,
                                                      final_model = "cubic",
                                                      periodicity = TRUE,
                                                      boot_trials = 100,
                                                      sim_trials = 100,
                                                      miss_bins = c(1,2,3,4,5),
                                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90)),
                           
                           
                           blasto_not_top2_baddley = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/blasto/",
                                                      base_path = "/Shared/Statepi_Diagnosis/prelim_results/blasto/",  # base path to original prelim extract results
                                                      out_path = "/Shared/Statepi_Diagnosis/projects/blasto/blasto_not_top2_baddley/",   # path to output delay new results
                                                      ssd_name = "blasto",
                                                      cp = 92,
                                                      upper_bound = 365,
                                                      final_model = "cubic",
                                                      periodicity = TRUE,
                                                      boot_trials = 100,
                                                      sim_trials = 100,
                                                      miss_bins = c(1,2,3,4,5),
                                                      duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90)),
                           
                           blasto_NA_state = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/blasto/",
                                                          base_path = "/Shared/Statepi_Diagnosis/prelim_results/blasto/",  # base path to original prelim extract results
                                                          out_path = "/Shared/Statepi_Diagnosis/projects/blasto/blasto_NA_state/",   # path to output delay new results
                                                          ssd_name = "blasto",
                                                          cp = 92,
                                                          upper_bound = 365,
                                                          final_model = "cubic",
                                                          periodicity = TRUE,
                                                          boot_trials = 100,
                                                          sim_trials = 100,
                                                          miss_bins = c(1,2,3,4,5),
                                                          duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60,90)),
                           
                           cocci = list(small_db_path = "/Shared/AML/truven_extracts/small_dbs/cocci/",
                                         base_path = "/Shared/Statepi_Diagnosis/prelim_results/cocci/",  # base path to original prelim extract results
                                         out_path = "/Shared/Statepi_Diagnosis/projects/cocci/",   # path to output delay new results
                                         ssd_name = "cocci",
                                         cp = 65,
                                         upper_bound = 365,
                                         final_model = "cubic",
                                         periodicity = TRUE,
                                         boot_trials = 100,
                                         sim_trials = 100,
                                         miss_bins = c(1,2,3,4,5),
                                         duration_bins = c(1,2,3,4,5,6,7,10,14,17,21,30,45,60))
)

# save parameters
save(final_delay_params,file = "/Volumes/AML/params/final_delay_params.RData")


