#!/bin/bash

source ./script-functions.sh

export PLIST=./model_lists/surf-performance
export ALIST=./model_lists/surf-analysis
export SDATABASE=surf.db
export SURF_DIR=surface_data

# steady state runs
if [ -z "$SKIP_STEADY_STATE" ]
then
    echo "Running steady state calcs..."
    skip_moles
    skip_analyt
    skip_flex
    skip_mass

    # get steady state times for all models
    declare -a ss_arr
    steady_args=-"-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
    ss_arr+=($(./launch.sh ./options/sp-opts single 1 1 $PLIST $steady_args  -S PlatinumLarge | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $ALIST $steady_args | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $ALIST  $steady_args --remove_falloff | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $ALIST  $steady_args --remove_thirdbody | grep -o -E '[0-9]{3,10}'))

    ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $ALIST  $steady_args --remove_falloff --remove_thirdbody | grep -o -E '[0-9]{3,10}'))

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
./launch.sh ./options/thresh-opts mpi 10 1 $ALIST -R performance -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $ALIST -R performance --remove_falloff -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $ALIST -R performance --remove_thirdbody -O $SURF_DIR -L
./launch.sh ./options/react-opts mpi 10 1 $ALIST -R performance --remove_thirdbody --remove_falloff -O $SURF_DIR -L
./launch.sh ./options/surf-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge

# analysis runs
./launch.sh ./options/thresh-opts single 1 1 $ALIST -R analysis -D $SDATABASE -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $ALIST -R analysis -D $SDATABASE --remove_falloff -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $ALIST -R analysis -D $SDATABASE --remove_thirdbody -O $SURF_DIR
./launch.sh ./options/react-opts single 1 1 $ALIST -R analysis -D $SDATABASE --remove_thirdbody --remove_falloff -O $SURF_DIR
./launch.sh ./options/surf-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
