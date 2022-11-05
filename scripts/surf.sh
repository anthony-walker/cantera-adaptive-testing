#!/bin/bash

#SBATCH --get-user-env                      #Use user env

#SBATCH -A nrg

#SBATCH -n 1                                # number of MPI tasks (default 1)

#SBATCH -c 20                               # cores

#SBATCH -p mime4							# name of partition or queue

#SBATCH --time=7-00:00:00

#SBATCH --mail-type=FAIL,END				# send email when job begins, ends or aborts

#SBATCH --mail-user=walkanth@oregonstate.edu		# send email to this address

#SBATCH -o surf-%j.out

# load any software environment module required for app

# run my jobs
echo "Slurm ID: $SLURM_JOB_ID"
echo "SURF_OPTION=$SURF_OPTION"
python surf-analysis.py $SURF_OPTION
