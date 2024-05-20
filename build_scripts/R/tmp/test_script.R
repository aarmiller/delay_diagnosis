
library(tidyverse)
library(smallDB)

args = commandArgs(trailingOnly=TRUE)

# name of condition
cond_name <- args[1]
# cond_name <- "histo"

con <- DBI::dbConnect(RSQLite::SQLite(), paste0("/Shared/AML/truven_extracts/small_dbs/",cond_name,"/",cond_name,".db"))

delay_base_path <- paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/")

load("/Shared/AML/params/delay_any_params.RData")
# load("/Volumes/AML/params/delay_any_params.RData")


out_path <- paste0(delay_base_path,"risk_models/")

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

delay_params <- delay_any_params[[cond_name]]
rm(delay_any_params)

paste0(out_path, "abx_duration_models.RData")


rmarkdown::render(input = "github/delay_diagnosis/build_scripts/R/report_scripts/abx_duration_report.Rmd",
                  params = list(cond = cond_name),
                  output_dir = out_path)