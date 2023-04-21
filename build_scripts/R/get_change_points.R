
args = commandArgs(trailingOnly=TRUE)

# devtools::install_github("aarmiller/smallDB", dependencies = FALSE, force = TRUE)
# devtools::install_github("aarmiller/codeBuildr", dependencies = FALSE, force = TRUE)
# codeBuildr::avail_ssd_codes()

# name of condition
cond_name <- args[1]
# cond_name <- "tb"

load("/Shared/AML/params/delay_any_params.RData")

delay_params <- delay_any_params[[cond_name]]

out_path <- paste0(delay_params$path,"change_point_results/")

if (!dir.exists(out_path)){
  dir.create(out_path)
}

rmarkdown::render(input = "github/delay_diagnosis/build_scripts/R/change_point_report.Rmd",
                  params = list(cond = cond_name),
                  output_dir = out_path)


