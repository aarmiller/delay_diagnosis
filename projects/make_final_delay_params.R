
################################################################################
### This Script Contains the parameters necessary for the final project analysis 
################################################################################

delay_any_params <- list(sarcoid = list(base_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/",
                                        cp = 150,
                                        upper_bound = 365*2,
                                        final_model = NA,
                                        add_validity_constraint = FALSE,
                                        miss_bins = c(1,2,3,4,5,7,10),
                                        duration_bins = c(3,7,14,21,30,60,90,120,150)),
                         
                         sarcoid_lung = list(base_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_lung/",
                                             cp = 150,
                                             upper_bound = 365*2,
                                             final_model = NA,
                                             add_validity_constraint = FALSE,
                                             miss_bins = c(1,2,3,4,5,7,10),
                                             duration_bins = c(3,7,14,21,30,60,90,120,150)),
                         
                         sarcoid_skin = list(base_path = "/Shared/Statepi_Diagnosis/projects/sarcoid/sarcoid_skin/",
                                             cp = 150,
                                             upper_bound = 365*2,
                                             final_model = NA,
                                             add_validity_constraint = FALSE,
                                             miss_bins = c(1,2,3,4,5,7,10),
                                             duration_bins = c(3,7,14,21,30,60,90,120,150))
                         
)

# save parameters
save(delay_any_params,file = "/Volumes/AML/params/delay_any_params.RData")


