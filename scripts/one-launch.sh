#!/bin/bash
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
source ./script-functions.sh
if [ -z "$1" ]
then
    echo "No arguements file given!"
    exit 1
else
    export ADD_ARGS="$(<$1)"
    DIRNAME=${ADD_ARGS##* }
    if [ ! -d $DIRNAME ]
    then
        mkdir $DIRNAME
    fi
fi
# current model to test
if [ -z "$2" ]
then
    echo "No model given!"
    exit 2
fi
# script file prefix
if [ -z "$3" ]
then
    echo "No script file prefix given!"
    exit 3
fi
# number of runs
if [ -z "$4" ]
then
    echo "No number of runs given!"
    exit 4
fi
# memory argument with a default
if [ -z "$5" ]
then
    export AMS="10G"
else
    export AMS=$5
fi
# loops in batch script
if [ -z "$6" ]
then
    echo "Number of batch loops not set, running 1 loop."
    export BATCH_LOOPS=1
else
    export BATCH_LOOPS=$4
fi
# make slurm-output dir
if [ ! -d "slurm-output" ]
then
    mkdir slurm-output
fi
# get runners
declare -a RUNNERS=()
define_runners "$3"
# run jobs
for job in "${RUNNERS[@]}"
do
    export CURR_MODEL=$2
    # evaluate jobs 10 times
    for (( i = 1; i <= ${4}; i++ ))
    do
        eval $job
    done
done

