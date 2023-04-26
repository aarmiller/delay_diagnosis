
args = commandArgs(trailingOnly=TRUE)

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
codeBuildr::avail_ssd_codes()

# name of condition
cond_name <- args[1]
# cond_name <- "tb"


rmarkdown::render(input = "github/delay_diagnosis/build_scripts/R/report_scripts/delay_report.Rmd",
                  params = list(cond = cond_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/"))

