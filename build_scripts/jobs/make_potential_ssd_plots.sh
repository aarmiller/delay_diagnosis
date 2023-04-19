#!/bin/bash

#Set the name of the job.
#$ -N make_potential_ssd_plots

#Set the shell that should be used to run the job.
#$ -S /bin/bash

#Set the directory for input an error files
#$ -o /localscratch/Users/aarmille/$JOB_NAME_o$JOB_ID.txt -e /localscratch/Users/aarmille/$JOB_NAME_e$JOB_ID.txt

#Select the queue to run in
#$ -q AML-HM

#Select the number of slots the job will use
#$ -pe smp 10


#Print information from the job into the output file
/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

#Description
desc="Make SSD plots for $1 "

#Script Paths
script_path="github/delay_diagnosis/build_scripts/R/make_potential_ssd_plots.R"
r_out="/Shared/AML/job_out/R_out/delay_jobs/make_potential_ssd_plots_$JOB_ID.txt"

# Print job info to job_history file
echo Job: $JOB_NAME "/" ID: $JOB_ID  "/" Date: `date` "/" Desc: $desc "/" Path: $script_path >> /Shared/AML/job_out/job_history/job_history.txt

#Send e-mail at beginning/end/suspension of job
#$ -m bes

#E-mail address to send to
#$ -M aaron-miller@uiowa.edu

#INPUT JOB HERE
module load stack/legacy
module load R
Rscript $script_path > $r_out $1

# Move the error and output files
mv $SGE_STDOUT_PATH /Shared/AML/job_out/SGE_out/delay_jobs
mv $SGE_STDERR_PATH /Shared/AML/job_out/SGE_out/delay_jobs
