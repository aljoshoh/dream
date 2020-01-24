#!/bin/sh
# #################################
# ##### JOB MASTER EXECUTABLE #####
# #################################

# TASKS:
#   Make and move directories
#   submit cluster job

# Define static variables
FILE=$3
tasksperarray=$6
BASEDIR=/storage/groups/cbm01/$4
BASEDIRarrayscript=/storage/groups/cbm01/tools/submit_job


# Check status of command
(($#==6)) || { echo -e "\nUSAGE: $0 <task-ID max> <task-ID-max> <file to execute> <working directory in cbm> <docker image> <tasks per array>\n\n"; exit; }
echo $BASEDIR
echo $FILE
filename="env_array_multi"
begintime="`date "+%d.%m.%Y - %H:%M:%S"`"


# Define non-static variables


#   # Nr. of parallel tasks
paralleltask1=$1
paralleltask2=$2


# Start Preprocessing
echo "Preprocessing..."
cd $BASEDIR #assure working directory is the directory of the file
echo Begintime: $begintime
echo "Finished preprocessing."
# End Preprocessing


############################
##### Launch Array Job #####
############################

echo "Launching job array..."
sbatch --output=/home/icb/$USER/slurm_output-%j.txt --error=/home/icb/$USER/slurm_error-%j.txt --ntasks=$6 --array=$paralleltask1-$paralleltask2 $BASEDIRarrayscript/$filename.sh "$BASEDIR" "$paralleltasks" "$FILE" "$4" "$5" "$6"

# TODO
# calculate paralelltasks because its not supplied
