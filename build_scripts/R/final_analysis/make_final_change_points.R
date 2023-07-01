
# This script is a placeholder for generating a change-point 

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
project_name <-  args[2]
cond_name <- "sarcoid"
project_name <-  "sarcoid_skin"
params = list(cond = cond_name,
              proj = project_name)

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_change_point_report.Rmd",
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/", cond_name, "/",
                                      project_name, "/"))
