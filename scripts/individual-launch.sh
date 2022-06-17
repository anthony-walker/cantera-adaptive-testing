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
for i in {0..99}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "IsoOctane" "ThreeMethylHeptane" "NHexadecane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
    do
        export CURR_MODEL=$j
        # Approximate preconditioner
        if [ -z "$SKIP_APPROX" ]
        then
            echo "Launching approx-$CURR_MODEL"
            sbatch -J "approx-$CURR_MODEL" ./slurm-batches/approx-precon-single.sh
            sleep 0.1
        fi
        # Analytical preconditioner
        if [ -z "$SKIP_ANALYT" ]
        then
        echo "Launching analyt-$CURR_MODEL"
        sbatch -J "analyt-$CURR_MODEL" ./slurm-batches/analyt-precon-single.sh
        sleep 0.1
        fi
    done
done

# Run mole and mass studies
for i in {0..99}
do
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    for j in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "IsoOctane" "ThreeMethylHeptane" "NHexadecane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
    do
        export CURR_MODEL=$j
        # Mass run
        if [ -z "$SKIP_MASS" ]
        then
            echo "Launching mass-$CURR_MODEL"
            sbatch -J "mass-$CURR_MODEL" ./slurm-batches/mass-single.sh
            sleep 0.1
        fi
        # Moles run
        if [ -z "$SKIP_MOLES" ]
        then
            echo "Launching moles-$CURR_MODEL"
            sbatch -J "moles-$CURR_MODEL" ./slurm-batches/moles-single.sh
            sleep 0.1
        fi
    done
done
