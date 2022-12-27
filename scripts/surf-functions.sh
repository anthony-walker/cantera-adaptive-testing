#!/bin/bash

source ./script-functions.sh

surface_performance_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=surf-performance
    # steady state runs
    if [ -z "$SKIP_STEADY_STATE" ]
    then
        echo "Running steady state calcs..."
        skip_moles
        skip_analyt
        skip_flex
        skip_mass
        # set start and end thresholds
        export TSTART=0
        export TEND=0
        # get steady state times for all models
        declare -a ss_arr
        steady_args=-"-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
        ss_arr+=($(./launch.sh ./options/sp-opts single 1 1 $PLIST $steady_args  -S PlatinumLarge | grep -o -E '[0-9]{3,10}'))
        # check if jobs are still active
        ss_running=true
        while [ $ss_running == true ]
        do
            ss_running=false
            OUTPUT=$(slurm_job_print)
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
    unset TEND
    unset TSTART
    skip_moles
    skip_analyt
    skip_flex
    # launch all runs now that ss is found.
    ./launch.sh ./options/sp-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge
}

surface_reaction_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=reaction.db
    export SURF_DIR=reaction-analysis
    # steady state runs
    if [ -z "$SKIP_STEADY_STATE" ]
    then
        echo "Running steady state calcs..."
        skip_moles
        skip_analyt
        skip_flex
        skip_mass
        # set start and end thresholds
        export TSTART=0
        export TEND=0
        # get steady state times for all models
        declare -a ss_arr
        steady_args=-"-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST $steady_args | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_falloff | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_thirdbody | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_falloff --remove_thirdbody | grep -o -E '[0-9]{3,10}'))
        # check if jobs are still active
        ss_running=true
        while [ $ss_running == true ]
        do
            ss_running=false
            OUTPUT=$(slurm_job_print)
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
    skip_flex
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_falloff -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_thirdbody -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_thirdbody --remove_falloff -O $SURF_DIR -L
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_falloff -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_thirdbody -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_thirdbody --remove_falloff -O $SURF_DIR
}

surface_threshold_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=threshold.db
    export SURF_DIR=threshold-analysis
    # steady state runs
    if [ -z "$SKIP_STEADY_STATE" ]
    then
        echo "Running steady state calcs..."
        skip_moles
        skip_analyt
        skip_flex
        skip_mass
        # set start and end thresholds
        export TSTART=0
        export TEND=0
        # get steady state times for all models
        declare -a ss_arr
        steady_args=-"-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST $steady_args | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_falloff | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_thirdbody | grep -o -E '[0-9]{3,10}'))
        ss_arr+=($(./launch.sh ./options/sa-opts single 1 1 $PLIST  $steady_args --remove_falloff --remove_thirdbody | grep -o -E '[0-9]{3,10}'))
        # check if jobs are still active
        ss_running=true
        while [ $ss_running == true ]
        do
            ss_running=false
            OUTPUT=$(slurm_job_print)
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
    skip_flex
    unset TSTART
    unset TEND
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
}
