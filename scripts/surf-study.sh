#!/bin/bash

source ./script-functions.sh

export MLIST=./model_lists/surf-models-copy
export SDATABASE=surf.db
export SURF_DIR=surface_data

# skip_moles
# skip_analyt
# skip_approx

# declare -a ss_arr

# ss_arr+=($(./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR | grep -o -E '[0-9]{3,10}'))

# ss_arr+=($(./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_falloff -O $SURF_DIR | grep -o -E '[0-9]{3,10}'))

# ss_arr+=($(./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_thirdbody -O $SURF_DIR | grep -o -E '[0-9]{3,10}'))

# ss_arr+=($(./launch.sh ./options/surf-opts single 1 1 $MLIST -R steady -MTS 1e-3 -MS 1e9 --remove_falloff --remove_thirdbody -O $SURF_DIR | grep -o -E '[0-9]{3,10}'))

# ss_running=true
# while [ $ss_running = true ]
# do
#     ss_running=false
#     OUTPUT=$(./job-print.sh)
#     for ss in "${ss_arr[@]}"
#     do
#         if [[ "$OUTPUT" == *"$ss"* ]]; then
#             echo "$ss: still found"
#             ss_running=true
#         fi
#     done
#     # sleep if ss is still runs
#     if [ $ss_running = true ]; then
#         echo "Steady state runs still active, sleeping for a minute..."
#         sleep 60
#     fi
# done

# skip certain runs
reset_skips
skip_moles
skip_analyt
skip_mass

# # performance runs
# ./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance -O $SURF_DIR -L
# ./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_falloff -O $SURF_DIR -L
# ./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_thirdbody -O $SURF_DIR -L
./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance --remove_thirdbody --remove_falloff -O $SURF_DIR -L

# # analysis runs
# ./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D $SDATABASE -O $SURF_DIR
# ./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D $SDATABASE --remove_falloff -O $SURF_DIR
# ./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D $SDATABASE --remove_thirdbody -O $SURF_DIR
# ./launch.sh ./options/surf-opts single 1 1 $MLIST -R analysis -D $SDATABASE --remove_thirdbody --remove_falloff -O $SURF_DIR
