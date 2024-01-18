
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
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh sepsis_revised10
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh sepsis_revised10

#### cvst ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh cvst
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh cvst
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh cvst

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh cvst
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh cvst


#### Sarcoid ####

# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh sarcoid
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh sarcoid
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh sarcoid


# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh sarcoid
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh sarcoid

#### Meningitis ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh meningitis
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh meningitis
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh meningitis

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh meningitis
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh meningitis
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh meningitis
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh meningitis
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh meningitis
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh meningitis
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh meningitis



#### Blasto ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh blasto
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh blasto
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh blasto

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh blasto



#### COCCI ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh cocci
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh cocci
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh cocci

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh cocci

#### PCP ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh pcp
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh pcp
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh pcp

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh pcp
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh pcp

#### Dengue ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh dengue
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh dengue
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh dengue

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh dengue


#### Pertussis ####
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh pertussis
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh pertussis
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh pertussis

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh pertussis
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh pertussis
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh pertussis


