

red_book <- codeBuildr::load_redbook()

red_book %>% glimpse()

red_book %>% count(THRDTDS, name = "Number of NDCs") %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/redbook_therapeutic_detailed_group.csv")

red_book %>% 
  count(GENNME, name = "Number of NDCs") %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/redbook_generic_name.csv")

red_book %>% 
  count(THRCLDS, name = "Number of NDCs") %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/redbook_therapeutic_class.csv")

red_book %>% 
  filter(str_detect(tolower(THRDTDS),"prednisone")) %>% glimpse()


red_book %>% 
  filter(str_detect(tolower(THRDTDS),"prednisone") |
           str_detect(tolower(PRODNME),"prednisone") | 
           str_detect(tolower(GENNME),"prednisone"))


find_codes <- function(term){
  red_book %>% 
    filter(str_detect(tolower(THRDTDS),term) |
             str_detect(tolower(PRODNME),term) | 
             str_detect(tolower(GENNME),term))
}

find_codes("prednisone")
find_codes("dexamethasone") %>% glimpse()

tmp <- tibble(term = c("prednisone","cortisone","dexamethasone","betamethasone","ethamethasoneb",
                "hydrocortisone","methylprednisolone","prednisolone","deflazacort",
                "triamcinolone","fludrocortisone","mometasone")) %>% 
  mutate(codes = map(term,find_codes))

tmp %>% 
  unnest(codes) %>% 
  select(-term) %>% 
  distinct() %>% 
  filter(ROADS=="Oral") %>% 
  select(NDCNUM,PKSIZE,GNINDDS:PRDCTDS,ORGBKDS:GENNME,ROADS) %>%  
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/oral_steroids.csv")

red_book %>% 
  filter(str_detect(tolower(THRDTDS),"term") |
                 +              str_detect(tolower(PRODNME),"term") | 
                 +              str_detect(tolower(GENNME),"term"))


library(codeBuildr)

codeBuildr::avail_ssd_codes("sarcoid")

codeBuildr::load_ssd_codes("sarcoid") %>% 
  rename(dx = code)  %>% 
  mutate(dx_ver = ifelse(type == "icd9",9L,10L)) %>% 
  select(-type) %>% 
  inner_join(codeBuildr::all_icd_labels) %>% 
  arrange(dx_ver,dx) %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/ssd_list.csv")


#### Inhaler Codes ####

inhaler_codes <- read_csv("~/OneDrive - University of Iowa/WorkingPapers/sepsis_delay_dx/data/all_inhaler_codes.csv")

inhaler_codes %>% glimpse()

inhaler_codes <- inhaler_codes %>% 
  select(NDCNUM,PKSIZE,GNINDDS:PRDCTDS,ORGBKDS:GENNME,ROADS)

inhaler_codes %>% 
  write_csv("~/OneDrive - University of Iowa/delay_dx_projects/sarcoid/code_sets/inhaler_codes.csv")



