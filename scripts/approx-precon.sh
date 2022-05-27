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

# starting threshold
if [ -z "$TSTART" ]
then
    TSTART=0
fi
# ending threshold
if [ -z "$TEND" ]
then
    TEND=18
fi

for th in $(seq $TSTART $TEND)
do
    if [ $th -ne 0 ]
    then
        THR="1e-$th"
    else
        THR="0"
    fi
    export PRECON_OPTS="$CURR_MODEL -L -v -M -P -T $THR $ADD_ARGS"
    echo $PRECON_OPTS
    mpirun -n 10 -hosts=$HOSTNAME adaptive-testing.mpi_run_same $PRECON_OPTS
    sleep 0.1
done
