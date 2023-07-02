
# This script is a placeholder for generating a change-point 

args = commandArgs(trailingOnly=TRUE)

# name of condition
project_name <- args[1]
# project_name <-  "sarcoid_lung"

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_change_point_report.Rmd",
                  params = list(proj = project_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/", 
                                      stringr::str_split(project_name, "_")[[1]][1], "/",
                                      project_name, "/"))
