

tmp1 <- read_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/data/inhaler_codes_combined/inhaler1.csv")
tmp2 <- read_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/data/inhaler_codes_combined/inhaler2.csv")
tmp3 <- read_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/data/inhaler_codes_combined/inhaler3.csv")

red_book <- codeBuildr::load_redbook()

bind_rows(tmp1 %>% select(NDCNUM),
          tmp2 %>% select(NDCNUM),
          tmp3 %>% select(NDCNUM)) %>% 
  distinct(NDCNUM) %>%
  mutate(NDCNUM = str_pad(NDCNUM,pad = "0",width = 11)) %>% 
  inner_join(red_book) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/data/all_inhaler_codes.csv")

red_book %>% 
  select(NDCNUM) %>% 
  mutate(len = nchar(NDCNUM)) %>% 
  count(len)
