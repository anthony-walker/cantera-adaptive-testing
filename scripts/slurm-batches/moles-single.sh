#!/bin/bash

#SBATCH --get-user-env                      #Use user env

#SBATCH -A nrg

#SBATCH -n 1                                # number of MPI tasks (default 1)

#SBATCH -p mime4, share								# name of partition or queue

#SBATCH --time=7-00:00:00

#SBATCH --mail-type=FAIL,END				# send email when job begins, ends or aborts

#SBATCH --mail-user=walkanth@oregonstate.edu		# send email to this address

#SBATCH -o ./slurm-output/%x-%j.out

# load any software environment module required for app

# run my jobs
# run my jobs
echo "Slurm ID: $SLURM_JOB_ID"
export MOLE_OPTS="$CURR_MODEL -L -v -M $ADD_ARGS"
echo "Options: $MOLE_OPTS"
# run mass
adaptive-testing $MOLE_OPTS
