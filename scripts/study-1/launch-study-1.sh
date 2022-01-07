#!/bin/bash
# use adaptive-test.sh to test launch
# update opts file for different options
mkdir slurm-output
mkdir data
for model in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "NHexadecane" "IsoOctane" "ThreeMethylHeptane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
do
    export CURR_MODEL=$model
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    export ADD_ARGS="$(<opts)"
    for i in {0..9}
    do
        sbatch -J "$model-prec-$i" precon.sh
        sleep0.1
        sbatch -J "$model-mass-$i" mass.sh
        sleep 0.1
	sbatch -J "$model-mole-$i" moles.sh
    done
done


