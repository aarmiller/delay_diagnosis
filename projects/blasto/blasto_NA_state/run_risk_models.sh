#!/bin/bash

#Set the name of the job.
#$ -N run_risk_models

#Set the shell that should be used to run the job.
#$ -S /bin/bash

#Set the directory for input an error files
#$ -o /localscratch/Users/atlan/$JOB_NAME_o$JOB_ID.txt -e /localscratch/Users/atlan/$JOB_NAME_e$JOB_ID.txt

#Select the queue to run in
#$ -q AML-HM

#Select the number of slots the job will use
#$ -pe smp 80

#Print information from the job into the output file
/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

#Description
desc="run_risk_models"

#Script Paths
script_path="/Shared/Statepi_Diagnosis/atlan/github/delay_diagnosis/projects/blasto/blasto_NA_state/run_risk_models.R"
r_out="/Shared/Statepi_Diagnosis/atlan/out/delay_diagnosis/run_risk_models_$JOB_ID.txt"

# Print job info to job_history file
echo Job: $JOB_NAME "/" ID: $JOB_ID  "/" Date: `date` "/" Desc: $desc "/" Path: $script_path >> /Shared/Statepi_Diagnosis/atlan/job_history/delay_diagnosis/job_history.txt

#Send e-mail at beginning/end/suspension of job
#$ -m bes

#E-mail address to send to
#$ -M alan-arakkal@uiowa.edu

#INPUT JOB HERE
module load stack/2020.2
module load r/4.0.2_gcc-8.4.0
Rscript $script_path > $r_out 

# Move the error and output files
mv $SGE_STDOUT_PATH /Shared/Statepi_Diagnosis/atlan/out/delay_diagnosis
mv $SGE_STDERR_PATH /Shared/Statepi_Diagnosis/atlan/out/delay_diagnosis