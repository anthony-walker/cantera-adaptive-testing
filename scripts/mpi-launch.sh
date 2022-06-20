#!/bin/bash
# get slurm script functions
source ./slurm-scripts.sh
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
check_args_and_dirs $1
# Run preconditioned studies because they are faster
declare -a RUNNERS=()
define_runners "mpi"
# run jobs
for job in "${RUNNERS[@]}"
do
    for i in {0..9}
    do
        # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
        for j in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "IsoOctane" "ThreeMethylHeptane" "NHexadecane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
        do
            export CURR_MODEL=$j
            eval $job
        done
    done
done


