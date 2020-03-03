#!/bin/bash

(($#==3)) || { echo -e "\nUSAGE: $0 <script name in scripts/ folder> <how many instances> <how many containers>\n\n"; exit; }

scriptname=$1
numberinstances=$2
container=$3

/storage/groups/cbm01/workspace/dream_aml/submit_job/execute_calc_slurm_multi_update.sh 1 $container $scriptname workspace/dream_aml/scripts r3.6.1_dream $numberinstances


