#!/bin/bash

#### A script for generating the cmd.list file. 
# Each line in this file will be an instance of a run_GutFunFind.py command
# to be run on each of the 4644 representative UHGG genomes under:
# /fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/*
#
# the only basic thing to do for this script is replace the function name 
# with your function of interest and this script will generate a list of commands
# to run in parallel on the cluster with script <>.sh

echo "- initializing variables"
DATA_BASE_DIR="data"
FUNCTION_NAME="Cysteine_Degradation"

echo "- grabbing each of the genome paths"
all_genomes=$(ls -d /fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/*/) # list of all genome folders
for genome_folder in $all_genomes
do
    genome_sub_folder=$(ls -d $genome_folder*/) # list of all genomes under the current genome folder
    for genome in $genome_sub_folder
    do
        genome_name=$(basename -- "$genome") # actual name of the current genome to be appended to the next path
        genome+="genome/$genome_name" # great! now $genome contains the path needed for the `-g` flag in GFF

        # echo "-- assembling the current line"
        echo "run_GutFunFind.py rep -b $DATA_BASE_DIR -f $FUNCTION_NAME -g $genome -o results/$FUNCTION_NAME/$genome_name/$genome_name" >> scripts/GFF_on_UGHH_rep_genomes.txt
    done
done

echo "- DONE"
