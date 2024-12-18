
rm(list=ls()[!(ls() %in% c("proj_name", "project_name", "cps_index", "cps_vector", "delay_params_org"))])

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
# codeBuildr::avail_ssd_codes()

# Set working directory
# name of condition
# args = commandArgs(trailingOnly=TRUE)
# project_name <- args[1]
# project_name <- "dengue"

load("/Shared/AML/params/final_delay_params.RData")
delay_params <- final_delay_params[[proj_name]]

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_delay_report.Rmd",
                  params = list(proj = project_name),
                  output_dir = delay_params$out_path)

