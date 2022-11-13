#!/bin/bash

source ./script-functions.sh

export MLIST=./model_lists/test-list

skip_moles
skip_analyt
skip_approx

./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9

reset_skips
skip_moles
skip_analyt

# ./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis

# ./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance
