#!/bin/bash

#Set the name of the job.
#$ -N make_delay_report

#Set the shell that should be used to run the job.
#$ -S /bin/bash

#Set the directory for input an error files
#$ -o /localscratch/Users/aarmille/$JOB_NAME_o$JOB_ID.txt -e /localscratch/Users/aarmille/$JOB_NAME_e$JOB_ID.txt

#Select the queue to run in
#$ -q AML-HM

#Select the number of slots the job will use
#$ -pe smp 30


#Print information from the job into the output file
/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

#Description
desc="Make delay base data for $1"

#Script Paths
script_path="github/truven_db_extracts/R/delay_scripts/make_delay_report.R"
r_out="/Shared/AML/job_out/R_out/db_extracts/make_delay_report_$JOB_ID.txt"

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
mv $SGE_STDOUT_PATH /Shared/AML/job_out/SGE_out/db_extracts
mv $SGE_STDERR_PATH /Shared/AML/job_out/SGE_out/db_extracts
