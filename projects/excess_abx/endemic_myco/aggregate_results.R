
rm(list = ls())

###################
#### Load Data ####
###################

#### Cocci ---------------------------------------------------------------------

load("~/Data/projects/excess_abx/cocci/data/report_data.RData")
load("~/Data/projects/excess_abx/cocci/data/baseline_data.RData")
cocci_res <- list(baseline_data = baseline_data,
                  weekly_cp = weekly_cp,
                  p1 = p1,
                  p2 = p2,
                  p3 = p3,
                  p4 = p4,
                  p5 = p5,
                  abx_name_weekly_fits = abx_name_weekly_fits,
                  cubic_res_plots = cubic_res_plots,
                  cubic_stats = cubic_stats,
                  lm_res_plots = lm_res_plots,
                  lm_stats = lm_stats,
                  quad_res_plots = quad_res_plots,
                  quad_stats = quad_stats,
                  excess_res1 = excess_res1,
                  excess_res2 = excess_res2,
                  index_dates = index_dates,
                  final_abx_counts =final_abx_counts,
                  agg_boot_res = agg_boot_res)

rm(weekly_cp,p1,p2,p3,p4,p5,abx_name_weekly_fits,cubic_res_plots,cubic_stats,
   lm_res_plots,lm_stats,quad_res_plots,quad_stats,excess_res1,excess_res2,
   index_dates,final_abx_counts,agg_boot_res,baseline_data)

#### Histo ---------------------------------------------------------------------

load("~/Data/projects/excess_abx/histo/data/report_data.RData")
load("~/Data/projects/excess_abx/histo/data/baseline_data.RData")
histo_res <- list(baseline_data = baseline_data,
                  weekly_cp = weekly_cp,
                  p1 = p1,
                  p2 = p2,
                  p3 = p3,
                  p4 = p4,
                  p5 = p5,
                  abx_name_weekly_fits = abx_name_weekly_fits,
                  cubic_res_plots = cubic_res_plots,
                  cubic_stats = cubic_stats,
                  lm_res_plots = lm_res_plots,
                  lm_stats = lm_stats,
                  quad_res_plots = quad_res_plots,
                  quad_stats = quad_stats,
                  excess_res1 = excess_res1,
                  excess_res2 = excess_res2,
                  index_dates = index_dates,
                  final_abx_counts =final_abx_counts,
                  agg_boot_res = agg_boot_res)

rm(weekly_cp,p1,p2,p3,p4,p5,abx_name_weekly_fits,cubic_res_plots,cubic_stats,
   lm_res_plots,lm_stats,quad_res_plots,quad_stats,excess_res1,excess_res2,
   index_dates,final_abx_counts,agg_boot_res,baseline_data)

#### Blasto --------------------------------------------------------------------

load("~/Data/projects/excess_abx/blasto/data/report_data.RData")
load("~/Data/projects/excess_abx/blasto/data/baseline_data.RData")
blasto_res <- list(baseline_data = baseline_data,
                   weekly_cp = weekly_cp,
                  p1 = p1,
                  p2 = p2,
                  p3 = p3,
                  p4 = p4,
                  p5 = p5,
                  abx_name_weekly_fits = abx_name_weekly_fits,
                  cubic_res_plots = cubic_res_plots,
                  cubic_stats = cubic_stats,
                  lm_res_plots = lm_res_plots,
                  lm_stats = lm_stats,
                  quad_res_plots = quad_res_plots,
                  quad_stats = quad_stats,
                  excess_res1 = excess_res1,
                  excess_res2 = excess_res2,
                  index_dates = index_dates,
                  final_abx_counts =final_abx_counts,
                  agg_boot_res = agg_boot_res)

rm(weekly_cp,p1,p2,p3,p4,p5,abx_name_weekly_fits,cubic_res_plots,cubic_stats,
   lm_res_plots,lm_stats,quad_res_plots,quad_stats,excess_res1,excess_res2,
   index_dates,final_abx_counts,agg_boot_res,baseline_data)


#### Sporotrichosis ------------------------------------------------------------

load("~/Data/projects/excess_abx/sporotrichosis/data/report_data.RData")
load("~/Data/projects/excess_abx/sporotrichosis/data/baseline_data.RData")
sporo_res <- list(baseline_data = baseline_data,
                  weekly_cp = weekly_cp,
                  p1 = p1,
                  p2 = p2,
                  p3 = p3,
                  p4 = p4,
                  p5 = p5,
                  abx_name_weekly_fits = abx_name_weekly_fits,
                  cubic_res_plots = cubic_res_plots,
                  cubic_stats = cubic_stats,
                  lm_res_plots = lm_res_plots,
                  lm_stats = lm_stats,
                  quad_res_plots = quad_res_plots,
                  quad_stats = quad_stats,
                  excess_res1 = excess_res1,
                  excess_res2 = excess_res2,
                  index_dates = index_dates,
                  final_abx_counts =final_abx_counts,
                  agg_boot_res = agg_boot_res)

rm(weekly_cp,p1,p2,p3,p4,p5,abx_name_weekly_fits,cubic_res_plots,cubic_stats,
   lm_res_plots,lm_stats,quad_res_plots,quad_stats,excess_res1,excess_res2,
   index_dates,final_abx_counts,agg_boot_res,baseline_data)

#########################
#### Combine Results ####
#########################


### Baseline stats -------------------------------------------------------------

rename(blasto_res$baseline_data, blasto = out) %>% 
  left_join(rename(cocci_res$baseline_data, cocci = out)) %>% 
  left_join(rename(histo_res$baseline_data, histo = out)) %>% 
  left_join(rename(sporo_res$baseline_data, sporo = out)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/baseline_stats.csv")



### Number and percent of antibiotics in excess --------------------------------

total_excess_abx <- select(blasto_res$excess_res1,
       abx, blasto = lm) %>% 
  full_join(select(cocci_res$excess_res1,
                    abx, cocci = lm)) %>% 
  full_join(select(histo_res$excess_res1,
                    abx, histo = lm) ) %>% 
  full_join(select(sporo_res$excess_res1,
                    abx, sporo = lm)) %>% 
  arrange(abx) %>% 
  filter(abx!="total") %>% 
  mutate_all(~replace_na(.,""))

total_excess_abx %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/total_excess_abx.csv")

pct_excess_abx <- select(blasto_res$excess_res2,
                         abx, blasto = lm) %>% 
  full_join(select(cocci_res$excess_res2,
                   abx, cocci = lm)) %>% 
  full_join(select(histo_res$excess_res2,
                   abx, histo = lm) ) %>% 
  full_join(select(sporo_res$excess_res2,
                   abx, sporo = lm)) %>% 
  arrange(abx) %>% 
  filter(abx!="total") %>% 
  mutate_all(~replace_na(.,""))

pct_excess_abx %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/pct_excess_abx.csv")

### Number of patients with excess antibiotic ----------------------------------

blasto_res$lm_stats$main_excess_stats %>%
  mutate_at(vars(mean:hi),~round(.,2)) %>% 
  mutate(blasto = paste0(mean, " (", lo,", ",hi,")")) %>% 
  select(measure,blasto) %>% 
  inner_join(cocci_res$lm_stats$main_excess_stats %>%
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(cocci = paste0(mean, " (", lo,", ",hi,")")) %>% 
               select(measure,cocci)) %>% 
  inner_join(histo_res$lm_stats$main_excess_stats %>%
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(histo = paste0(mean, " (", lo,", ",hi,")")) %>% 
               select(measure,histo)) %>% 
  inner_join(sporo_res$lm_stats$main_excess_stats %>%
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(sporo = paste0(mean, " (", lo,", ",hi,")")) %>% 
               select(measure,sporo)) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/main_patient_stats.csv")


blasto_res$lm_stats$main_excess_stats_by_type %>% 
  filter(measure == "Percent of Patients with Excess ABX") %>% 
  mutate_at(vars(mean:hi),~round(.,2)) %>% 
  mutate(blasto = paste0(mean,"\n (",lo,", ",hi,")")) %>% 
  select(name,blasto) %>% 
  full_join(cocci_res$lm_stats$main_excess_stats_by_type %>% 
               filter(measure == "Percent of Patients with Excess ABX") %>% 
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(cocci = paste0(mean,"\n (",lo,", ",hi,")")) %>% 
               select(name,cocci)) %>% 
  full_join(histo_res$lm_stats$main_excess_stats_by_type %>% 
               filter(measure == "Percent of Patients with Excess ABX") %>% 
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(histo = paste0(mean,"\n (",lo,", ",hi,")")) %>% 
               select(name,histo)) %>% 
  full_join(sporo_res$lm_stats$main_excess_stats_by_type %>% 
               filter(measure == "Percent of Patients with Excess ABX") %>% 
               mutate_at(vars(mean:hi),~round(.,2)) %>% 
               mutate(sporo = paste0(mean,"\n (",lo,", ",hi,")")) %>% 
               select(name,sporo)) %>% 
  mutate_all(~replace_na(.,"")) %>% 
  arrange(name) %>% 
  write_csv("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/pct_w_excess_by_type.csv")


###############
#### Plots ####
###############

## Antibiotic trends by type ---------------------------------------------------

# Blasto
blasto_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = mean_fit), color = "red", size = 1) +
  geom_ribbon(aes(ymin = lo_fit, ymax = hi_fit), fill = "red", alpha = 0.2) +
  geom_line(aes(y = lo_fit), color = "red", linetype = 2) +
  geom_line(aes(y = hi_fit), color = "red", linetype = 2) +
  geom_vline(aes(xintercept = blasto_res$weekly_cp), linetype = 2) +
  facet_wrap(~name, scale = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Before Index Blastomycosis Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/blasto_expected.pdf")

blasto_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>%
  filter(week<=15) %>% 
  ggplot(aes(week,excess)) +
  geom_point() +
  scale_x_reverse() +
  geom_pointrange(aes(ymin = excess_low, ymax = excess_high)) +
  facet_wrap(~name, scale = "free") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Number of Excess Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/blasto_excess.pdf")

# Cocci
cocci_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = mean_fit), color = "red", size = 1) +
  geom_ribbon(aes(ymin = lo_fit, ymax = hi_fit), fill = "red", alpha = 0.2) +
  geom_line(aes(y = lo_fit), color = "red", linetype = 2) +
  geom_line(aes(y = hi_fit), color = "red", linetype = 2) +
  geom_vline(aes(xintercept = cocci_res$weekly_cp), linetype = 2) +
  facet_wrap(~name, scale = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Before Index Coccidioidomycosis Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/cocci_expected.pdf")

cocci_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>%
  filter(week<=15) %>% 
  ggplot(aes(week,excess)) +
  geom_point() +
  scale_x_reverse() +
  geom_pointrange(aes(ymin = excess_low, ymax = excess_high)) +
  facet_wrap(~name, scale = "free") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Number of Excess Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/cocci_excess.pdf")

# Histo
histo_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = mean_fit), color = "red", size = 1) +
  geom_ribbon(aes(ymin = lo_fit, ymax = hi_fit), fill = "red", alpha = 0.2) +
  geom_line(aes(y = lo_fit), color = "red", linetype = 2) +
  geom_line(aes(y = hi_fit), color = "red", linetype = 2) +
  geom_vline(aes(xintercept = histo_res$weekly_cp), linetype = 2) +
  facet_wrap(~name, scale = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Before Index Histoplasmosis Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/histo_expected.pdf")

histo_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>%
  filter(week<=15) %>% 
  ggplot(aes(week,excess)) +
  geom_point() +
  scale_x_reverse() +
  geom_pointrange(aes(ymin = excess_low, ymax = excess_high)) +
  facet_wrap(~name, scale = "free") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Number of Excess Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/histo_excess.pdf")

# Sporo
sporo_res$agg_boot_res %>% 
  filter(model == "lm") %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = mean_fit), color = "red", size = 1) +
  geom_ribbon(aes(ymin = lo_fit, ymax = hi_fit), fill = "red", alpha = 0.2) +
  geom_line(aes(y = lo_fit), color = "red", linetype = 2) +
  geom_line(aes(y = hi_fit), color = "red", linetype = 2) +
  geom_vline(aes(xintercept = sporo_res$weekly_cp), linetype = 2) +
  facet_wrap(~name, scale = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Before Index Sporotrichosis Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/sporo_expected.pdf")

sporo_res$agg_boot_res %>% 
  filter(!(name %in% c("total","total (individual sum)"))) %>%
  filter(week<=15) %>% 
  ggplot(aes(week,excess)) +
  geom_point() +
  scale_x_reverse() +
  geom_pointrange(aes(ymin = excess_low, ymax = excess_high)) +
  facet_wrap(~name, scale = "free") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylab("Number of Excess Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/sporo_excess.pdf")

## Antibiotic trends total -----------------------------------------------------

bind_rows(blasto_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Coccidioidomycosis",
                   cp = blasto_res$weekly_cp),
          cocci_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Blastomycosis",
                   cp = cocci_res$weekly_cp),
          histo_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Histoplasmosis",
                   cp = histo_res$weekly_cp),
          sporo_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Sporotrichosis",
                   cp = sporo_res$weekly_cp)) %>% 
  ggplot(aes(week,n)) +
  geom_point() +
  scale_x_reverse() +
  geom_line(aes(y = mean_fit), color = "red", size = 1) +
  geom_ribbon(aes(ymin = lo_fit, ymax = hi_fit), fill = "red", alpha = 0.2) +
  geom_line(aes(y = lo_fit), color = "red", linetype = 2) +
  geom_line(aes(y = hi_fit), color = "red", linetype = 2) +
  geom_vline(aes(xintercept = cp), linetype = 2) +
  facet_wrap(~cond, scale = "free_y") +
  theme_bw() +
  ylab("Number of Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/total_expected.pdf")

bind_rows(blasto_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Coccidioidomycosis",
                   cp = blasto_res$weekly_cp),
          cocci_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Blastomycosis",
                   cp = cocci_res$weekly_cp),
          histo_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Histoplasmosis",
                   cp = histo_res$weekly_cp),
          sporo_res$agg_boot_res %>% 
            filter(model == "lm") %>% 
            filter(name=="total (individual sum)") %>% 
            mutate(cond = "Sporotrichosis",
                   cp = sporo_res$weekly_cp)) %>% 
  filter(week<=cp) %>% 
  ggplot(aes(week,excess)) +
  geom_point() +
  scale_x_reverse() +
  geom_pointrange(aes(ymin = excess_low, ymax = excess_high)) +
  facet_wrap(~cond, scale = "free") +
  theme_bw() +
  ylab("Number of Excess Antibiotic Prescriptions") +
  xlab("Weeks Before Index Diagnosis")
ggsave("~/OneDrive - University of Iowa/WorkingPapers/excess_abx/endemic_myco/results/total_excess.pdf")

