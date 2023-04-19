

library(tidyverse)

load("/Volumes/AML/params/delay_any_params.RData")



tmp <- codeBuildr::avail_disease_codes() %>% 
  inner_join(tibble(name = names(delay_any_params))) %>% 
  mutate(cp = map_dbl(name,~delay_any_params[[.]]$cp)) %>% 
  mutate(upper_bound = map_dbl(name,~delay_any_params[[.]]$upper_bound)) %>% 
  mutate(final_model = map_chr(name,~delay_any_params[[.]]$final_model)) 


for (i in 1:nrow(tmp)){
  if (!dir.exists(paste0("projects/",tmp$name[i]))) {
    dir.create(paste0("projects/",tmp$name[i]))
  }
  
  if (!file.exists(paste0("projects/",tmp$name[i],paste0("/readme.md")))){
    
    fileConn<-file(paste0("projects/",tmp$name[i],paste0("/readme.md")))
    writeLines(c(paste0("# Project sub-directory for ",tmp$description[i]),
                 paste0("This sub-directory contains scripts for conducting the final analysis for ", tmp$description[i]),
                 "",
                 "## Scripts",
                 "The following is a list and summary of scripts contained in this directory:",
                 "",
                 "",
                 "## Sub-projects",
                 paste0("The following is a list and description of any sub-projects conducted for ",tmp$description[i],":")), fileConn)
    close(fileConn) 
  }
}

for (i in 1:nrow(tmp)){
  if (!dir.exists(paste0("results/projects/",tmp$name[i]))) {
    dir.create(paste0("results/projects/",tmp$name[i]))
  }
  
  if (!file.exists(paste0("results/projects/",tmp$name[i],paste0("/readme.md")))){
    
    fileConn<-file(paste0("results/projects/",tmp$name[i],paste0("/readme.md")))
    writeLines(c(paste0("# Results sub-directory for ",tmp$description[i]),
                 paste0("This sub-directory contains results for the final analysis of ", tmp$description[i]),
                 ""), fileConn)
    close(fileConn) 
  }
}
