library(tidyverse)

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/functions/make_abx_plots_fun.R")

conds <- list.dirs("/Shared/Statepi_Diagnosis/prelim_results", recursive = F, full.names = F)

tmp <- parallel::mclapply(conds, 
                          function(x){make_abx_plots(cond_name = x)},
                          mc.cores = 30)



make_abx_plots(cond_name = "sporotrichosis")
