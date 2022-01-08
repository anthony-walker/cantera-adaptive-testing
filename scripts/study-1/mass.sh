#!/bin/bash

#SBATCH --get-user-env                      #Use user env

#SBATCH -A nrg

#SBATCH -n 10                                # number of MPI tasks (default 1)

#SBATCH -p mime4								# name of partition or queue

#SBATCH --time=7-00:00:00

#SBATCH --mail-type=FAIL,END				# send email when job begins, ends or aborts

#SBATCH --mail-user=walkanth@oregonstate.edu		# send email to this address

#SBATCH -o ./slurm-output/%x-%j.out

# load any software environment module required for app

# run my jobs
echo "Slurm ID: $SLURM_JOB_ID"
export MASS_OPTS="$CURR_MODEL -L -v $ADD_ARGS"
echo "Options: $MASS_OPTS"
# run mass
mpirun -n 10 -hosts=$HOSTNAME adaptive-testing.mpi_run_same $MASS_OPTS