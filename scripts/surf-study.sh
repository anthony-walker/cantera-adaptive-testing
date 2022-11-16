#!/bin/bash

source ./script-functions.sh

export MLIST=./model_lists/surf-models

skip_moles
skip_analyt
skip_approx

./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9

./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_falloff

./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_thirdbody

./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_falloff --remove_thirdbody

reset_skips
skip_moles
skip_analyt

# normal runs
./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D "surf-data.db"

./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance

# remove falloff runs
./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D "surf-data.db" --remove_falloff

./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_falloff

# remove thirdbody runs
./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D "surf-data.db" --remove_thirdbody

./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_thirdbody

# remove both runs
./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D "surf-data.db" --remove_thirdbody --remove_falloff

./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_thirdbody --remove_falloff
