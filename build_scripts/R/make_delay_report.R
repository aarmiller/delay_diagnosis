
args = commandArgs(trailingOnly=TRUE)

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
codeBuildr::avail_ssd_codes()

# name of condition
cond_name <- args[1]
# cond_name <- "tb"

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

rmarkdown::render(input = "github/truven_db_extracts/R/delay_scripts/delay_report.Rmd",
                  params = list(cond = cond_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/"))

