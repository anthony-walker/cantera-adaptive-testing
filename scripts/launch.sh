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
# Run preconditioned studies because they are faster
for i in {0..9}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","4G" "Butane","4G" "NHeptane","4G" "IsoOctane","4G" "ThreeMethylHeptane","6G" "NHexadecane","6G" "MethylFiveDeconate","6G" "MethylDeconateNHeptane","8G" "TwoMethylnonadecane","12G"
    do
        IFS=","; 
        set -- $j;
        export model=$1
        export AMS=$2
        sbatch -J "$model-prec-$i" precon.sh --mem=$AMS
        sleep 0.1
    done
done

# Run mole and mass studies
for i in {0..9}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","4G" "Butane","4G" "NHeptane","4G" "IsoOctane","4G" "ThreeMethylHeptane","6G" "NHexadecane","6G" "MethylFiveDeconate","6G" "MethylDeconateNHeptane","8G" "TwoMethylnonadecane","16G"
    do
        IFS=","; 
        set -- $j;
        export model=$1
        export AMS=$2
        sbatch -J "$model-mass-$i" mass.sh --mem=$AMS
        sleep 0.1
	    sbatch -J "$model-mole-$i" moles.sh --mem=$AMS
        sleep 0.1
    done
done