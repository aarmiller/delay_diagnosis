library(tidyverse)

all_cond <- tibble(cond = c("endocarditis", "epidural_abscess", "meningitis", "sepsis", "vertebral_osteo"),
                    file_cond = c("endoca", "epidur", "mening", "sepsis", "verteb") )

# cond <- "sepsis"

for (conds in all_cond$cond){
  
 
  out_path <- paste0("/Shared/Statepi_Diagnosis/projects/", conds, "_kaiser/temperature_exploration/")
  
  if(!dir.exists(out_path)){
    dir.create(out_path, recursive = T)
  }
  
  tryCatch({
    
  rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/kaiser/temperature_exploration/temperature_time_of_day_blood_culture.Rmd",
                    output_dir = out_path,
                    params = list(cond = conds))
  
  }, error = function(e) {
    
    message("An error occurred while rendering the Rmd for : ", conds)
    
  })
  
  rm(list = ls()[!ls() %in%  c("all_cond", "conds")])
  
}