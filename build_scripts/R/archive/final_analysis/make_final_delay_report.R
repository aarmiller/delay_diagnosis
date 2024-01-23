
args = commandArgs(trailingOnly=TRUE)

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
# codeBuildr::avail_ssd_codes()

# Set working directory
# name of condition
project_name <- args[1]

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/final_analysis/final_delay_report.Rmd",
                  params = list(proj = project_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/", 
                                      stringr::str_split(project_name, "_")[[1]][1], 
                                      ifelse( stringr::str_split(project_name, "_")[[1]][1] == project_name ,
                                              "/", paste0("/", project_name, "/")))
                  )

