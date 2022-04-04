#!/bin/bash
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
if [ -z "$1" ]
then
    export ADD_ARGS=""
    mkdir data
else
    export ADD_ARGS="$(<$1)"
    mkdir ${ADD_ARGS##* }
fi
# make slurm-output dir
mkdir slurm-output

# Run mole and mass studies
for i in {0..1}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "MethylDeconateNHeptane","20G"
    do
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        sbatch -J "$CURR_MODEL-mass-$i" mass.sh --mem=$AMS
        sleep 0.1
	sbatch -J "$CURR_MODEL-mole-$i" moles.sh --mem=$AMS
        sleep 0.1
    done
done
