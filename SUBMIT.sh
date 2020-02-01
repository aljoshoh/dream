#!/bin/bash

(($#==2)) || { echo -e "\nUSAGE: $0 <script name in scripts/ folder> <how many instances>\n\n"; exit; }

sciptname=$1
numberinstances=$2

/storage/groups/cbm01/workspace/dream_aml/submit_job/execute_calc_slurm.sh 1 1 $scriptname storage/groups/cbm01/workspace/dream_aml/scipts dream_aml $numberinstances

