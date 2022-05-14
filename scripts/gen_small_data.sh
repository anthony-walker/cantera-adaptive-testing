#!/bin/bash
for i in {1..2}
do
    for mod in "Hydrogen" "MethaneGRI"
    do
        # preconditioned
        adaptive-testing $mod -L -M -P -T 1e-8
        adaptive-testing $mod -L -M -P -T 1e-12
        adaptive-testing $mod -L -M -P -T 0
        adaptive-testing $mod -L -M -P -T 1e-9
        # mole reactor no preconditioning
        adaptive-testing $mod -L -M
        # no moles
        adaptive-testing $mod -L
    done
done



