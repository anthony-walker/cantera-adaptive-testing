#!/bin/bash

#SBATCH --get-user-env                      #Use user env

#SBATCH -A nrg

#SBATCH -n 10                                # number of MPI tasks (default 1)

#SBATCH -p mime4,share								# name of partition or queue

#SBATCH --time=7-00:00:00

#SBATCH --mail-type=FAIL,END				# send email when job begins, ends or aborts

#SBATCH --mail-user=walkanth@oregonstate.edu		# send email to this address

#SBATCH -o ./slurm-output/%x-%j.out

# load any software environment module required for app

# run my jobs
echo "Slurm ID: $SLURM_JOB_ID"
# zero threshold
export PRECON_OPTS="$CURR_MODEL -L -v -M -P -S GMRES -T 0 $ADD_ARGS"
echo $PRECON_OPTS
mpirun -n 10 -hosts=$HOSTNAME adaptive-testing.mpi_run_same $PRECON_OPTS
sleep 0.1
# varying thresholds
for th in {1..18}
do
    echo $PRECON_OPTS
    export PRECON_OPTS="$CURR_MODEL -L -v -M -P -S GMRES -T 1e-$th $ADD_ARGS"
    mpirun -n 10 -hosts=$HOSTNAME adaptive-testing.mpi_run_same $PRECON_OPTS
    sleep 0.1
done
