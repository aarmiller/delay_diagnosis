library(tidyverse)

rm(list=ls()[!(ls() %in% c("cond_name"))])
condition_name <- stringr::str_to_title(str_remove(cond_name, "_kaiser"))

rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/delay_report.Rmd",
                  params = list(cond = cond_name,
                                condition_name = condition_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/prelim_results/",cond_name,"/delay_results/"))

