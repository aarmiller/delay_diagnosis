
# This script is a placeholder for generating a change-point 


# name of condition
project_name <- "dengue_validated"


rmarkdown::render(input = "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/dengue/dengue_validated/final_change_point_report.Rmd",
                  params = list(proj = project_name),
                  output_dir = paste0("/Shared/Statepi_Diagnosis/projects/", 
                                      stringr::str_split(project_name, "_")[[1]][1], 
                                      ifelse( stringr::str_split(project_name, "_")[[1]][1] == project_name ,
                                              "/", paste0("/", project_name, "/")))
                  )
