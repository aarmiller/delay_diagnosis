

cond <- "sepsis"
out_path <- paste0("/Shared/Statepi_Diagnosis/projects/", cond, "_kaiser/temperature_exploration/")

if(!dir.exists(out_path)){
  dir.create(out_path, recursive = T)
}


rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/kaiser/temperature_exploration/temperature_blood_culture.Rmd",
                  output_dir = out_path,
                  params = list(cond = cond))