#!/bin/bash
# get slurm script functions
source ./script-functions.sh
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
check_args_and_dirs $1
# mpi or single arg
if [ -z "$2" ]
then
    echo "No run type specified, use mpi or single"
    exit 2
else
    export RTYPE=$2
fi
# number of loops
# mpi or single arg
if [ -z "$3" ]
then
    echo "Number of loops not set, running 1 loop."
    export LOOPS=1
else
    export LOOPS=$3
fi
# loops in batch script
if [ -z "$4" ]
then
    echo "Number of batch loops not set, running 1 loop."
    export BATCH_LOOPS=1
else
    export BATCH_LOOPS=$4
fi
# Run preconditioned studies because they are faster
declare -a RUNNERS=()
declare -a MODELS=(`cat models`)
define_runners "$RTYPE"
# run jobs
for job in "${RUNNERS[@]}"
do
    for j in "${MODELS[@]}"
    do
        # break model arguments up
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        # evaluate jobs 10 times
        for (( i = 1; i <= ${LOOPS}; i++ ))
        do
            eval $job
        done
    done
done


