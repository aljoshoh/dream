#!/bin/bash

#SBATCH -J arrayscript
#SBATCH --mem-per-cpu=3G
#SBATCH -p icb_cpu
#SBATCH -t 2-00:00:00
#SBATCH --cpus-per-task 1

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

# this function is deprecated
create_dir () {
	directory="$1"
	root="$2"
	if [ -d "$root/$directory" ]; then
		directory=$(($directory+1))
		directory=$(create_dir $directory $root)
	fi
	mkdir -p /localscratch/$USER/$directory
	echo Created directory in lscratch: $USER/$directory
}

# update the number of tasks
#array=$SLURM_JOB_NODELIST
#newtasks=${#array[@]}
newtasks=$SLURM_JOB_NUM_NODES
echo Number of new tasks: $SLURM_JOB_NUM_NODES
numtasks=$(($newtasks+$tasksperarray+$newtasks))
echo Ntasks set to $numtasks !!!
scontrol update JobId=$SLURM_JOB_ID NumTasks=$numtasks #add tasks for adding and removing images 


# execute these lines for each node we are working on...
echo Nodelist: $SLURM_JOB_NODELIST
echo Directory in lscratch: $SLURM_JOB_ID
dir=$SLURM_JOB_ID
for node in $SLURM_JOB_NODELIST;
do
	echo Executing create-docker on $node
	srun --nodelist=$node -n $newtasks --ntasks-per-node=1 bash -c 'DOCKER=$0; dir=$SLURM_JOB_ID; mkdir -p /localscratch/$USER/$SLURM_JOB_ID;
	cd /localscratch/$USER/$dir;
	echo Created directory for local environment: $dir, with docker image: $DOCKER	
ch-tar2dir /storage/groups/cbm01/tools/alexander.ohnmacht/r-studio-charliecloud-master-b76b94e8c5040cd58fb0ffd1a463fd7409bb886e/exports/$DOCKER.tar.gz /localscratch/$USER/$dir;
	echo Set up the docker image !;
	' $DOCKER
done
wait

##############################
### The actual task submission
##############################

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
	srun -n 1 ch-run -b /localscratch/:/localscratch/ -b /storage/groups/:/storage/groups/ /localscratch/$USER/$dir/$DOCKER/ -- Rscript $BASEDIR/$FILE $i $TASKID &
	sleep 5
done

wait

echo delete charlieclould image
for node in $SLURM_JOB_NODELIST;
do
        echo Executing delete-docker on $node
	srun --nodelist=$node -n $newtasks --ntasks-per-node=1 bash -c 'rm -rf /localscratch/$USER/$dir/'
done
wait

# i am done
echo The job has ended.

