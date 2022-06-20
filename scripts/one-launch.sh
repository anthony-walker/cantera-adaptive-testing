#!/bin/bash
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
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
    AMS="10G"
else
    AMS=$5
fi
# make slurm-output dir
if [ ! -d "slurm-output" ]
then
    mkdir slurm-output
fi
# run over number of runs
for j in $(seq 1 $4)
do
    SCRIPT="./slurm-batches/$3-single.sh"
    export CURR_MODEL=$2
    sbatch -J "$2-$3" $SCRIPT --mem=$AMS
done

