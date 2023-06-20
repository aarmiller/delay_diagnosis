
#### tb ####

### main extracts (Note, these are run before the delay diagnosis scripts, this generates the database on the AML drive)
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh tb
# qsub github/truven_db_extracts/jobs/main_scripts/build_index_db_enroll_restrict_valid.sh tb 365
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_enroll_restrict_valid.sh tb 365

### delay_diagnosis scripts
# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots_restricted.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh tb

#### histo ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh histo
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh histo
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh histo
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh histo

### delay_diagnosis scripts
# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh tb


#### Sepsis ####
# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh sepsis_revised10
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh sepsis_revised10
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh sepsis_revised10
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh sepsis_revised10