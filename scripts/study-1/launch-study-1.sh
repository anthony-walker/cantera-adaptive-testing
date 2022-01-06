#!/bin/bash
# use adaptive-test.sh to test launch
# update opts file for different options
mkdir data   
for model in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "NHexadecane" "IsoOctane" "ThreeMethylHeptane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
do  
    export CURR_MODEL=$model
    # ADD_ARGS is a string argumenent given to the the launch file to add arguements to the run e.g. "--no_vol_prob --no_net_prob" or "--max_time_step 1e-8"
    export ADD_ARGS="$(<opts)" 
    for i in {0..9}
    do
        sbatch -J "adaptive-prec-$model-$i" adaptive-precon.sh
        sbatch -J "adaptive-mass-$model-$i" adaptive-mass.sh
        sbatch -J "adaptive-mole-$model-$i" adaptive-moles.sh
    done
done


