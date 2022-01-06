#!/bin/bash
for i in {1..2}
do
    for mech in "Hydrogen" "MethaneGRI"
    do
        # preconditioned
        adaptive-testing $mech -L -M -P -S GMRES -T 1e-8 --no_net_prob
        # mole reactor no preconditioning
        adaptive-testing $mech -L -M --no_net_prob
        # no moles
        adaptive-testing $mech -L --no_net_prob
    done
done



