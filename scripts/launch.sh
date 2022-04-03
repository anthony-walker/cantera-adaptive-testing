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
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","4G" "Butane","4G" "NHeptane","6G" "IsoOctane","8G" "ThreeMethylHeptane","12G" "NHexadecane","12G" "MethylFiveDeconate","16G" "MethylDeconateNHeptane","20G" "TwoMethylnonadecane","24G"
    do
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        echo "Launching $CURR_MODEL-prec-approx-$i -- $AMS"
        sbatch -J "$CURR_MODEL-prec-$i" approx-precon.sh --mem=$AMS
        sleep 0.1
        echo "Launching $CURR_MODEL-prec-analyt-$i -- $AMS"
        sbatch -J "$CURR_MODEL-prec-$i" analyt-precon.sh --mem=$AMS
        sleep 0.1
    done
done

# Run mole and mass studies
for i in {0..9}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","4G" "Butane","4G" "NHeptane","6G" "IsoOctane","8G" "ThreeMethylHeptane","12G" "NHexadecane","12G" "MethylFiveDeconate","16G" "MethylDeconateNHeptane","20G" "TwoMethylnonadecane","24G"
    do
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        echo "Launching $CURR_MODEL-mass-$i -- $AMS"
        sbatch -J "$CURR_MODEL-mass-$i" mass.sh --mem=$AMS
        sleep 0.1
        echo "Launching $CURR_MODEL-moles-$i -- $AMS"
	    sbatch -J "$CURR_MODEL-mole-$i" moles.sh --mem=$AMS
        sleep 0.1
    done
done
