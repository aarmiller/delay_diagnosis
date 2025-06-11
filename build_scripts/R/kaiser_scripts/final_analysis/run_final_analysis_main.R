
# Run all final_analysis scripts

args = commandArgs(trailingOnly=TRUE)
project_name <- args[1]
proj_name <- project_name

load("/Shared/AML/params/final_delay_params_kaiser.RData")

delay_params_org <- final_delay_params
cps_vector <- final_delay_params[[project_name]]$cp

for (cps_index in cps_vector){
  
  print(paste("cps =", cps_index, "started"))
  
  final_delay_params <- delay_params_org
  final_delay_params[[project_name]]$cp <- cps_index
  final_delay_params[[project_name]]$out_path <- paste0(final_delay_params[[project_name]]$out_path, "delay_window_1_", cps_index - 1, "/")
  
  if (!dir.exists(final_delay_params[[project_name]]$out_path)) {
    dir.create(final_delay_params[[project_name]]$out_path)
  }
  
  save(final_delay_params, file = "/Shared/AML/params/final_delay_params_kaiser.RData")
  
  ## Run final_sim
  source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/final_sim.R")

  ## Run final_analyze_sim_res
  source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/final_analyze_sim_res.R")

  ## Run make_final_delay_report
  source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/make_final_delay_report.R")
  
  print(paste("cps =", cps_index, "finished"))
  
}

final_delay_params <- delay_params_org
save(final_delay_params,file = "/Shared/AML/params/final_delay_params_kaiser.RData")

## Run make_potential_ssd_plots
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/make_potential_ssd_plots.R")

## Run make_final_change_points
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/make_final_change_points.R")
