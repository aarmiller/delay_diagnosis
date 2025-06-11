

args = commandArgs(trailingOnly=TRUE)
cond_name <- args[1]

# run make_delay_base_data.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/make_delay_base_data.R")

# run make_potential_ssd_plots.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/make_potential_ssd_plots.R")

# run get_change_points.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/get_change_points.R")

# run get_delay_res_any.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/get_delay_res_any.R")

# run get_delay_res_ssd.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/get_delay_res_ssd.R")

# run make_delay_report.R ---------------------------------------------------

source("/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/R/kaiser_scripts/make_delay_report.R")