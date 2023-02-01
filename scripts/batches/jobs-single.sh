#!/bin/bash

#SBATCH --get-user-env                      #Use user env

#SBATCH -A nrg

#SBATCH -n 1                                # number of MPI tasks (default 1)

#SBATCH -p mime4,share							# name of partition or queue

#SBATCH --time=7-00:00:00

#SBATCH --mail-type=FAIL				# send email when job begins, ends or aborts

#SBATCH --mail-user=walkanth@oregonstate.edu		# send email to this address

#SBATCH -o ./slurm-output/%x-%j.out

#SBATCH --constraint=haswell

#SBATCH --exclusive

# load any software environment module required for app

# run my jobs
echo "Slurm ID: $SLURM_JOB_ID"
echo
echo $JOB_OPTIONS
echo
for i in {1..$BATCH_LOOPS}
do
    adaptive-testing $JOB_OPTIONS
    sleep $SLEEP_TIMER
done
