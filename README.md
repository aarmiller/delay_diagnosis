# Delayed Diagnosis Repository (delay_diagnosis)

This repository contains scripts and results for the various delayed diagnosis projects conducted by\... There are three primary folders associated with the various delayed diagnosis projects:

-   **build_scripts** - contains R build scripts and batch job scripts for performing the initial build and preliminary delayed diagnosis analysis for each disease
-   **projects** - contains the final analysis scripts for individual projects after a p
-   **results** - contains summary results for the various projects based on the preliminary analysis performed using the build scripts and final analysis performed using scripts in the projects folder
-   **manuscript material** - contains scripts for generating results (plots, tables, etc.) for each manuscript produced

## Build Scripts

The build scripts folder contains two subfolders: `R` and `jobs`. The `R` folder contains R scripts for performing various tasks and the `jobs` folder contains the batch scripts for performing these tasks. (Note, these build processes are performed once a small database extract has been created using extract scripts contained in the `truven_db_extracts` repository.) The following tasks are performed by scripts in this folder.

-   **make_potential_ssd_plots** - generates SSD curves for each of the top 100 ICD-9/10 codes up to the index diagnosis

-   **make_delay_base_data** - generates the required datasets for performing clustering and delayed diagnosis analysis

-   **get_clusters** - performs cluster analysis on different SSD trends before or after index diagnosis

-   **get_change_points** - performs change point analysis for

-   **get_delay_res_ssd** - runs a preliminary delay diagnosis analysis for SSD visits

-   **get_delay_res_any** - runs a preliminary delay diagnosis analysis for any visits

-   **make_delay_report** - generates a report of delayed diagnosis results

-   **run_risk_models** - runs a risk model analysis using a set of standard conditions

NOTE: Each of these scripts requires an entry into the `delay_any_params.RData` file located at `/AML/params/delay_any_params.RData`. This file contains the necessary parameters (e.g., data location, change-point, upper bound window, final model to use, etc.) for performing analysis. This file should be regularly updated throughout the analysis process. The script for generating and updating this file is located at `/delay_diagnosis/build_scripts/make_delay_params.R`. A corresponding file for performing the final analysis for each disease is also located within the projects folder.

## Build Scripts

A *project* defines the analysis and enrollment criteria corresponding to a particular study cohort. A project is created for each main condition (e.g., tuberculosis). Whenever a sub-population of a particular *main condition* is identified (e.g., pulmonary tuberculosis), a new project is created; such sub-projects are denoted by the naming convention *main condiion_sub-population* (e.g., tb_pulmonary).
