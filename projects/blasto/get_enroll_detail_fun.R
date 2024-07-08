collect_table <- function(settings = c("inpatient","outpatient","rx"),sources = c("ccae","mdcr"),
                          years=1:21, medicaid_years = NULL) {
  
  # convert years to stings
  years <- stringr::str_pad(years,2,pad="0")
  
  out <- tibble::tibble(setting=settings) %>%
    dplyr::mutate(source=purrr::map(.data$setting,~sources)) %>%
    tidyr::unnest(cols = c(source)) %>%
    dplyr::mutate(year=purrr::map(source,~years)) %>%
    tidyr::unnest(cols = c(year))
  
  if (!is.null(medicaid_years)) {
    medicaid_rows <- tibble::tibble(setting=settings) %>%
      dplyr::mutate(source="medicaid") %>%
      tidyr::unnest(cols = c(source)) %>%
      dplyr::mutate(year=purrr::map(source,~medicaid_years)) %>%
      tidyr::unnest(cols = c(year))
    
    out <- rbind(out,medicaid_rows)
  }
  
  return(out)
}


get_collapse_enrollment <- function (table, enrolid_list, collect_n = Inf,
                                     vars = c("egeoloc", "msa", "plantyp" ,"indstry"), db_path) {
  db_con <- src_sqlite(db_path)
  
  temp <- dplyr::tbl(db_con, paste0("enrollment_detail_", table$source, "_",  table$year)) %>%
    dplyr::filter(patient_id %in% enrolid_list) %>% 
    dplyr::select(c("patient_id", "dtstart", "dtend", vars)) %>% 
    dplyr::collect(n = collect_n) %>% 
    dplyr::mutate(patient_id = as.integer(.data$patient_id))
  
  out <- temp %>% dplyr::mutate_at(.vars = dplyr::vars(vars), .funs = list(as.integer))
  
  return(out)
}

gather_collapse_enrollment <- function (collect_tab = collect_table(), enrolid_list, collect_n = Inf,
                                        vars = c("egeoloc", "msa", "plantyp" ,"indstry"), 
                                        db_path, num_cores = NULL) {
  # require some pacakges
  require(tidyverse)
  require(dplyr)
  
  db_path2 <- db_path
  vars2 <- vars
  collect_n2 <- collect_n
  enrolid_list2 <- enrolid_list
  
  temp <- collect_tab %>% dplyr::select(-.data$setting)
  
  # set up clusters
  if (is.null(num_cores)) {
    num_cores <- min(nrow(temp), parallel::detectCores() - 1)
  } else {
    num_cores <- num_cores
  }
  
  cluster <- parallel::makeCluster(num_cores)
  parallel::clusterExport(cluster, varlist = c("get_collapse_enrollment"))
  
  parallel::clusterCall(cluster, function() library(tidyverse))
  parallel::clusterCall(cluster, function() library(dplyr))
  
  tmp <- parallel::parLapply(cl = cluster,
                             1:nrow(temp),
                             function(x){get_collapse_enrollment(table = temp[x, ],
                                                                 enrolid_list = enrolid_list2,
                                                                 vars = vars2,
                                                                 db_path = db_path2,
                                                                 collect_n = collect_n2)})
  
  parallel::stopCluster(cluster)
  gc()
  
  out <- tibble()
  for (i in 1:length(tmp)){
    x <- tmp[[i]] %>% nest(data = everything())
    out <- bind_rows(out, x)
  }
  temp <- out %>% select(data) %>% unnest() 
  
  temp_strata <- temp  %>% dplyr::select(c("patient_id", vars)) %>% 
    dplyr::distinct() %>% dplyr::mutate(strata = dplyr::row_number())
  
  temp <- temp %>% dplyr::inner_join(temp_strata, by = c("patient_id", vars))
  
  out <- temp %>% dplyr::arrange(.data$patient_id, .data$dtstart) %>% 
    dplyr::group_by(.data$patient_id) %>% 
    dplyr::mutate(gap =((.data$dtstart - dplyr::lag(.data$dtend)) > 1) | .data$strata != dplyr::lag(.data$strata), 
                  gap = ifelse(is.na(.data$gap), FALSE, .data$gap)) %>% 
    dplyr::mutate(period = cumsum(.data$gap)) %>% 
    dplyr::group_by_at(c("patient_id", "period", vars)) %>%
    dplyr::summarise(dtstart = min(.data$dtstart), 
                     dtend = max(.data$dtend)) %>% 
    dplyr::ungroup()
  return(out)
}