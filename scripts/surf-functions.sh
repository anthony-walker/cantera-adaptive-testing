#!/bin/bash

source ./script-functions.sh

surface_performance_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=surf-performance
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    unset TSTART
    unset TEND
    # launch all runs now that ss is found.
    ./launch.sh ./options/sp-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -E 1 -S PlatinumLarge
}

surface_reaction_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=reaction.db
    export SURF_DIR=reaction-analysis
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -E 1
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_falloff -O $SURF_DIR -L -E 1
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_thirdbody -O $SURF_DIR -L -E 1
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --remove_thirdbody --remove_falloff -O $SURF_DIR -L -E 1
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -E 1
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_falloff -O $SURF_DIR -E 1
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_thirdbody -O $SURF_DIR -E 1
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --remove_thirdbody --remove_falloff -O $SURF_DIR -E 1
}

surface_threshold_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=threshold.db
    export SURF_DIR=threshold-analysis
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    unset TSTART
    unset TEND
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -E 1
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -E 1
}
