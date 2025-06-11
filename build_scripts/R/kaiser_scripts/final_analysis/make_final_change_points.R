
# This script is a placeholder for generating a change-point 

rm(list=ls()[!(ls() %in% c("proj_name", "project_name"))])
# args = commandArgs(trailingOnly=TRUE)

# name of condition
# project_name <- args[1]
# project_name <-  "dengue"

load("/Shared/AML/params/final_delay_params_kaiser.RData")
delay_params <- final_delay_params[[proj_name]]

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/final_analysis/final_change_point_report.Rmd",
                  params = list(proj = project_name),
                  output_dir = delay_params$out_path)