#!/bin/bash

source ./script-functions.sh

surface_performance_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=performance_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=20
    # launch all runs now that ss is found.
    ./launch.sh ./options/sp-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge
}

surface_reaction_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=reaction.db
    export SURF_DIR=reaction_data
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L
    skip_mass
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_falloff -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody --enable_falloff -O $SURF_DIR -L
    # analysis runs
    unset SKIP_MASS
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
    skip_mass
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_falloff -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody --enable_falloff -O $SURF_DIR
}

jet_reaction_study() {
    # analysis runs
    export PLIST=./model_lists/jet-fuels
    export SDATABASE=jr.db
    export SURF_DIR=jr_data
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L
    skip_mass
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_falloff -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody -O $SURF_DIR -L
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody --enable_falloff -O $SURF_DIR -L
    # analysis runs
    unset SKIP_MASS
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
    skip_mass
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_falloff -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody -O $SURF_DIR
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody --enable_falloff -O $SURF_DIR
}

jet_reaction_fixed_timestep() {
    # analysis runs
    export PLIST=./model_lists/jet-fuels
    export SDATABASE=jrf.db
    export SURF_DIR=jrf_data
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -RS 1000
    skip_mass
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_falloff -O $SURF_DIR -L -RS 1000
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody -O $SURF_DIR -L -RS 1000
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody --enable_falloff -O $SURF_DIR -L -RS 1000
    # analysis runs
    unset SKIP_MASS
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -RS 1000
    skip_mass
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_falloff -O $SURF_DIR -RS 1000
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody -O $SURF_DIR -RS 1000
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody --enable_falloff -O $SURF_DIR -RS 1000
}

jet_threshold_study() {
    # analysis runs
    export PLIST=./model_lists/jet-fuels
    export SDATABASE=jth.db
    export SURF_DIR=jth_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=24
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L --enable_thirdbody --enable_falloff
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR --enable_thirdbody --enable_falloff
}


jet_threshold_fixed() {
    # analysis runs
    export PLIST=./model_lists/jet-fuels
    export SDATABASE=jft.db
    export SURF_DIR=jft_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=24
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -RS 1000 --enable_thirdbody --enable_falloff
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -RS 1000 --enable_thirdbody --enable_falloff
}

surface_threshold_study() {
    # analysis runs
    export PLIST=./model_lists/surf-analysis
    export SDATABASE=threshold.db
    export SURF_DIR=threshold_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=24
    # performance runs
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR
}

wsr_performance_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=performance_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=20
    # launch all runs now that ss is found.
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge
}

ttest() {
    factor=10
    for i in {1..20}
    do
        factor=$(echo "scale=$i; $factor*0.1" | bc)
        echo $factor
    done
}

run_all_surface_studies() {
    surface_threshold_study
    surface_reaction_study
    surface_performance_study
}
