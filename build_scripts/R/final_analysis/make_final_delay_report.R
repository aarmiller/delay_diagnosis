
args = commandArgs(trailingOnly=TRUE)

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
# codeBuildr::avail_ssd_codes()

# Set working directory
# name of condition
cond_name <- args[1]
project_name <-  args[2]
cond_name <- "sarcoid"
project_name <-  "sarcoid_skin"
params = list(cond = cond_name,
              proj = project_name)

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_delay_report.Rmd",
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/", cond_name, "/",
                                      project_name, "/"))

