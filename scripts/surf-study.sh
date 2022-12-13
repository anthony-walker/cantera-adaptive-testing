#!/bin/bash

source ./script-functions.sh

export MLIST=./model_lists/surf-models
export RLIST=./models_lists/reaction-analysis
export TLIST=./models_lists/thresh-analysis
export SDATABASE=surf.db
export SURF_DIR=surface_data

# steady state runs
if [ -z "$SKIP_STEADY_STATE" ]
then
    echo "Running steady state calcs..."
    skip_moles
    skip_analyt
    skip_flex
    skip_approx

    declare -a ss_arr
    steady_args=-"-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
    ss_arr+=($(./launch.sh ./options/surf-opts single 1 1 $MLIST $steady_args | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/thresh-opts single 1 1 $TLIST $steady_args | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/react-opts single 1 1 $RLIST  $steady_args --remove_falloff | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/react-opts single 1 1 $RLIST  $steady_args --remove_thirdbody | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/react-opts single 1 1 $RLIST  $steady_args --remove_falloff --remove_thirdbody | grep -o -E '[0-9]{3,10}'))

    ss_running=true
    while [ $ss_running == true ]
    do
        ss_running=false
        OUTPUT=$(./job-print.sh)
        for ss in "${ss_arr[@]}"
        do
            if [[ "$OUTPUT" == *"$ss"* ]]; then
                echo "$ss: still found"
                ss_running=true
            fi
        done
        # sleep if ss is still runs
        if [ $ss_running == true ]; then
            echo "Steady state runs still active, sleeping for a minute..."
            sleep 60
        fi
    done
fi

# skip certain runs
reset_skips
skip_moles
skip_analyt

# performance runs
./launch.sh ./options/surf-opts mpi 10 1 $MLIST -R performance -O $SURF_DIR -L
./launch.sh ./options/thresh-opts mpi 10 1 $TLIST -R performance -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $RLIST -R performance --remove_falloff -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $RLIST -R performance --remove_thirdbody -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $RLIST -R performance --remove_thirdbody --remove_falloff -O $SURF_DIR -L

# analysis runs
./launch.sh ./options/thresh-opts single 1 1 $TLIST -R analysis -D $SDATABASE -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $RLIST -R analysis -D $SDATABASE --remove_falloff -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $MLIST -R analysis -D $SDATABASE --remove_thirdbody -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $MLIST -R analysis -D $SDATABASE --remove_thirdbody --remove_falloff -O $SURF_DIR
