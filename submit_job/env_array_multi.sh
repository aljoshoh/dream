#!/bin/bash

#SBATCH -J arrayscript
#SBATCH --mem-per-cpu=5G
#SBATCH -p icb_cpu
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task 1
#SBATCH -N 1

BASEDIR=$1
TOTALTASKS=$2
FILE=$3
DOCKER=$5
tasksperarray=$6
TASKID=$SLURM_ARRAY_TASK_ID
SLEEPMULTIPLIER=10
#$5 is the docker image
echo $BASEDIR
echo Starting arrayscript job ...
echo USER: $USER

create_dir () {
	directory="$1"
	root="$2"
	if [ -d "$root/$directory" ]; then
		directory=$(($directory+1))
		directory=$(create_dir $directory $root)
	fi
	mkdir -p /localscratch/$USER/$directory
	echo $directory
}

dir=$(create_dir $TASKID /localscratch/$USER)
cd /localscratch/$USER/$dir

echo Created directory for local environment: $dir
ch-tar2dir /storage/groups/cbm01/tools/alexander.ohnmacht/r-studio-charliecloud-master-b76b94e8c5040cd58fb0ffd1a463fd7409bb886e/exports/$DOCKER.tar.gz /localscratch/$USER/$dir
echo Set up the docker image !


# some initial info
echo It is now:
date
echo
echo Running on machine
hostname
echo
echo Operating system
uname -r

# give system some rest to get the correct order
sleep $(( TASKID * SLEEPMULTIPLIER ))
echo "Launching task Nr. $TASKID out of $TOTALTASKS !"
sleep $(( TASKID * SLEEPMULTIPLIER ))
# sleep $(( TOTALTASKS * SLEEPMULTIPLIER ))

# code to execute with passing $TASKID to the R script
for i in `seq 1 $tasksperarray`;
do
	echo Executing $BASEDIR/$FILE with option $i
	srun -n 1 ch-run -b /storage/groups/:/storage/groups/ /localscratch/$USER/$dir/$DOCKER/ -- Rscript $BASEDIR/$FILE $i $TASKID &
	sleep 5
done

wait
echo delete charlieclould image
srun rm -rf /localscratch/$USER/$dir/$DOCKER/

# i am done
echo The job has ended.

