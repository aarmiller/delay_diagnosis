
# Run all final_analysis scripts

args = commandArgs(trailingOnly=TRUE)
project_name <- args[1]
proj_name <- project_name

## Run final_sim
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_sim.R")

## Run final_analyze_sim_res
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_analyze_sim_res.R")

## Run make_potential_ssd_plots
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/make_potential_ssd_plots.R")

## Run make_final_change_points
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/make_final_change_points.R")

## Run make_final_delay_report
source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/make_final_delay_report.R")