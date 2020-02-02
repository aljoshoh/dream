#!/bin/bash

(($#==2)) || { echo -e "\nUSAGE: $0 <script name in scripts/ folder> <how many instances>\n\n"; exit; }

scriptname=$1
numberinstances=$2

/storage/groups/cbm01/workspace/dream_aml/submit_job/execute_calc_slurm_multi.sh 1 1 $scriptname workspace/dream_aml/scripts r3.6.1_dream $numberinstances

