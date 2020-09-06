#!/bin/bash
#SBATCH -J run_GFF_cmds # Job name (short for bilirubun pipeline)
#SBATCH -o logs_2/run_GFF_cmds_%j_%a.o # output file
#SBATCH -e logs_2/run_GFF_cmds_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=08:00:00
#SBATCH --qos=throughput
#SBATCH --mem=16gb
#SBATCH --array=1-4644

##### [description]

# MAKE SURE THIS COMMAND GETS RUN IN A TERMINAL BEFORE RUNNING THIS SCRIPT
#bash ./scripts/generate_cmd_list.sh

# ALSO- MAKE SURE THAT YOUR results/FUNCTION_NAME/ folder is empty before running
#rm -rf results/Cysteine_Degradation/*
#rm -rf logs/*

# adding blast+ module
module purge && module load blast+
wait

# get current command from list of commands
current_cmd=`head -n ${SLURM_ARRAY_TASK_ID} scripts/GFF_on_UGHH_rep_genomes.txt | tail -n 1`

# make folder in results/ for each representative genome
echo $current_cmd | rev | cut -d '/' -f 1 | rev
mkdir "results/Cysteine_Degradation_2/$(echo $current_cmd | rev | cut -d '/' -f 1 | rev)"

# running list of commands
eval $current_cmd
