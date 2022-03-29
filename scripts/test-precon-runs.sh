#!/bin/bash

# Run preconditioned studies because they are faster

for j in "Hydrogen" "MethaneGRI" "DME" "JetA" "Butane" "NHeptane" "IsoOctane" "ThreeMethylHeptane" "NHexadecane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
do
    echo $j
    adaptive-testing $j -O test-prec-data -P -L -M -S GMRES -v
done


