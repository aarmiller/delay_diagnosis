
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
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh dengue
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh dengue
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh dengue

qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh dengue
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


#### endocarditis ####
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx_new.sh endocarditis
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh endocarditis
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh endocarditis


qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh endocarditis



#### blasto ####
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh blasto
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh blasto
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts_new.sh blasto


qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_clusters.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh endocarditis
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh endocarditis



#### Dengue ####
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh dengue
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh dengue
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh dengue;

qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh dengue       # Need to check this with new algorithm for delay days
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh dengue 
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh dengue
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh dengue


qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh cftr_main_cohort;

# Run these on 5/1/24

# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh dengue      
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh dengue 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh dengue
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh dengue

# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh chf  # Check this is done

# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh sarcoid
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh sarcoid

# Run these on 5/2/24

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh sarcoid       
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh sarcoid 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh sarcoid

# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh cocci
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh cocci
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh cocci

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh cocci       
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh cocci 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh cocci

# Run these on 5/3/24

# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh cocci

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh chf

# Run on 5/9/24


# Needed to fix index dates
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh cocci
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh cocci

# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh histo
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh histo

# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh blasto

######################
## Run on 5/14/2024 ##
######################

# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh blasto

# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh ami


# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh cocci       
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh cocci 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh cocci

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh histo       
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh histo 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh histo

# qsub github/delay_diagnosis/build_scripts/jobs/make_potential_ssd_plots.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh blasto       
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh blasto 
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh blasto


# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh ami

# fix # delay stuff AMI
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh ami# delay stuff AMI
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh ami
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh sarcoid
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh sarcoid

# fix lung cancer
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh lung_cancer
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh lung_cancer

# fix pe
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh pe
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh pe

########################
### Run on 5/15/2024 ###
########################

# fix chf
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh chf
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh chf

# delay stuff AMI
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh ami
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh ami
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh ami
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh ami
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh ami
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh ami

# delay stuff sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh sarcoid

# delay stuff lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh lung_cancer
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh lung_cancer

# delay stuff pe
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh pe
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh pe
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh pe
# qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh pe
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh pe
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh pe

# pull ANCA Vasculitis


#### Run on 5/20/24 ####

# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh histo
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh blasto
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh cocci
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh ami
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh sarcoid
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh pe
# qsub github/delay_diagnosis/build_scripts/jobs/run_duration_models.sh lung_cancer


# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh tb
# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh tb
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh tb
# qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh tb
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_any.sh tb
qsub github/delay_diagnosis/build_scripts/jobs/get_delay_res_ssd.sh tb
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_report.sh tb
qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh tb


# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh pertussis
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh endocarditis
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh append
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh hsv_enceph
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh nontb_myco
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh cmv
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh ebv
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh dvt
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh interstitial_lung
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh non_hodgkins_lymphoma
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh bladder_cancer
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh venous_insuf
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh ulcerative_colitis
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh crohns
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh hodgkins_lymphoma
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh renal_cell_cancer
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh acute_myeloid_leukemia
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh hairy_cell_leukemia

# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh kawasaki
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh sle
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh wegeners

# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh ra
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh takayasu
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh periodic_fever
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh adult_stills
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh thyroiditis
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh behcet
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh giant_cell_arteritis
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh epidural_abs
# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh sepsis

# qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh sinusitis_chronic

# qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh cftr_sinusitis_chronic
# qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh cftr_sinusitis_chronic

## Run on 6/11/24

# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh lung_cancer
# qsub github/delay_diagnosis/build_scripts/jobs/run_risk_models.sh pe

# Run on 9/16/24

qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh blasto
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh cocci
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh sarcoid



qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh hiv
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh bronchiectasis
qsub github/truven_db_extracts/jobs/main_scripts/get_index_dx.sh pcp


qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh hiv
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh bronchiectasis
qsub github/truven_db_extracts/jobs/main_scripts/build_small_db.sh pcp

qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh hiv
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh bronchiectasis
qsub github/truven_db_extracts/jobs/main_scripts/get_all_visit_counts.sh pcp


qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh hiv
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh bronchiectasis
qsub github/delay_diagnosis/build_scripts/jobs/make_delay_base_data.sh pcp


qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh blasto
qsub github/delay_diagnosis/build_scripts/jobs/get_change_points.sh cocci
