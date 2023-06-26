
# Inclusion Procedures
valid_procs <- c("31628",
                 "31624",
                 "31629",
                 "3324",
                 "31622",
                 "31620",
                 "31623",
                 "31633",
                 "31632",
                 "3323",
                 "3322",
                 "31645",
                 "31625",
                 "88172",
                 "3324",
                 "10022",
                 "76942",
                 "32405",
                 "38510",
                 "38505",
                 "3328",
                 "38525",
                 "3326",
                 "4011")


#Add in proc validation
proc_test <- all_procs %>% filter(proc %in% valid_procs)

#See if test is within 30 days before index date
proc_enrollees <- proc_test %>%
  filter(days_since_index>=-30 & days_since_index<=0)