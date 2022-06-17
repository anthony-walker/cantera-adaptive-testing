#!/bin/bash
# use launch slurm scripts for testing
# pass opts as argument and make output data dir
if [ -z "$1" ]
then
    echo "No args file given"
    exit 1
else
    export ADD_ARGS="$(<$1)"
    DIRNAME=${ADD_ARGS##* }
    if [ ! -d $DIRNAME ]
    then
        mkdir $DIRNAME
    fi
fi
# make slurm-output dir
if [ ! -d "slurm-output" ]
then
    mkdir slurm-output
fi
# Run preconditioned studies because they are faster
for i in {0..9}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","2G" "Butane","2G" "NHeptane","4G" "IsoOctane","4G" "ThreeMethylHeptane","4G" "NHexadecane","6G" "MethylFiveDeconate","8G" "MethylDeconateNHeptane","10G" "TwoMethylnonadecane","10G"
    do
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        # Approximate preconditioner
        if [ -z "$SKIP_APPROX" ]
        then
            echo "Launching approx-$CURR_MODEL -- $AMS"
            sbatch -J "approx-$CURR_MODEL" ./slurm-batches/approx-precon.sh --mem=$AMS
            sleep 0.1
        fi
        # Analytical preconditioner
        if [ -z "$SKIP_ANALYT" ]
        then
            echo "Launching analyt-$CURR_MODEL -- $AMS"
            sbatch -J "analyt-$CURR_MODEL" ./slurm-batches/analyt-precon.sh --mem=$AMS
            sleep 0.1
        fi
    done
done
# Run mole and mass studies
for i in {0..9}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen","2G" "MethaneGRI","2G" "DME","2G" "JetA","2G" "Butane","2G" "NHeptane","4G" "IsoOctane","4G" "ThreeMethylHeptane","4G" "NHexadecane","6G" "MethylFiveDeconate","8G" "MethylDeconateNHeptane","10G" "TwoMethylnonadecane","10G"
    do
        IFS=",";
        set -- $j;
        export CURR_MODEL=$1
        export AMS=$2
        # Mass run
        if [ -z "$SKIP_MASS" ]
        then
            echo "Launching mass-$CURR_MODEL -- $AMS"
            sbatch -J "mass-$CURR_MODEL" ./slurm-batches/mass.sh --mem=$AMS
            sleep 0.1
        fi
        # Moles run
        if [ -z "$SKIP_MOLES" ]
        then
            echo "Launching moles-$CURR_MODEL -- $AMS"
            sbatch -J "moles-$CURR_MODEL" ./slurm-batches/moles.sh --mem=$AMS
            sleep 0.1
        fi
    done
done
