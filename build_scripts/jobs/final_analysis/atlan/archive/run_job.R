
qsub "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/jobs/final_analysis/atlan/run_final_sim.sh" "dengue"

qsub "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/jobs/final_analysis/atlan/run_final_analyze_sim_res.sh" "dengue"

qsub "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/jobs/final_analysis/atlan/run_make_final_change_points.sh" "dengue"

qsub "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/jobs/final_analysis/atlan/run_make_potential_ssd_plots.sh" "dengue"

qsub "/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/build_scripts/jobs/final_analysis/atlan/run_make_final_delay_report.sh" "dengue"
