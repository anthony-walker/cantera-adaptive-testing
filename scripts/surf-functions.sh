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
    export TEND=18
    # launch all runs now that ss is found.
    ./launch.sh ./options/sp-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge -E 0.001
}

surface_short_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=perf_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    export TSTART=0
    export TEND=0
    # launch all runs now that ss is found.
    ./launch.sh ./options/nce-nab-pfr mpi 1 1 $PLIST -R performance -O $SURF_DIR -L -S PlatinumLarge -E 0.001 -ASF
}

perf_analysis_study() {
    # performance runs
    export PLIST=./model_lists/surf-performance
    export SURF_DIR=perf_data
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    skip_mass
    export TSTART=0
    export TEND=0
    # launch all runs now that ss is found.
    ./launch.sh ./options/pfr-opts single 1 1 $PLIST -R analysis -O $SURF_DIR -L -S PlatinumLarge -E 0.001 -D perf.db
}

network_series_test() {
    # analysis runs
    export PLIST=./model_lists/network-mods
    export SURF_DIR=series_data
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    # performance runs
    for i in {1..25}
    do
        ./launch.sh ./options/nrs-opts mpi 1 1 $PLIST -R performance -O $SURF_DIR -L -E 0.001 --nrs $i
        ./launch.sh ./options/nrs-opts mpi 1 1 $PLIST -R performance -O $SURF_DIR -L -E 0.001 --nrs $i --series
    done
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
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L -E 0.005
    skip_mass
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_falloff -O $SURF_DIR -L -E 0.005
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody -O $SURF_DIR -L -E 0.005
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody --enable_falloff -O $SURF_DIR -L -E 0.005
    # analysis runs
    unset SKIP_MASS
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -E 0.005
    skip_mass
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_falloff -O $SURF_DIR -E 0.005
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody -O $SURF_DIR -E 0.005
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody --enable_falloff -O $SURF_DIR -E 0.005
}

surface_condition_reaction_study() {
    # analysis runs
    export PLIST=./model_lists/jet-fuels
    export SDATABASE=jet_cond.db
    export SURF_DIR=jet_data
    export TSTART=0
    export TEND=0
    # skip certain runs
    reset_skips
    skip_moles
    skip_analyt
    skip_flex
    skip_mass
    # launch all runs now that ss is found.
    ./launch.sh ./options/pfr-opts single 1 1 $PLIST -R analysis -O $SURF_DIR -L -E 0.005 -D $SDATABASE
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
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L  -E 0.005
    skip_mass
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_falloff -O $SURF_DIR -L -E 0.005
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody -O $SURF_DIR -L -E 0.005
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance --enable_thirdbody --enable_falloff -O $SURF_DIR -L -E 0.005
    # analysis runs
    unset SKIP_MASS
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR -E 0.005
    skip_mass
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_falloff -O $SURF_DIR -E 0.005
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody -O $SURF_DIR -E 0.005
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE --enable_thirdbody --enable_falloff -O $SURF_DIR -E 0.005
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
    ./launch.sh ./options/sa-opts mpi 10 1 $PLIST -R performance -O $SURF_DIR -L --enable_thirdbody --enable_falloff -E 0.005
    # analysis runs
    ./launch.sh ./options/sa-opts single 1 1 $PLIST -R analysis -D $SDATABASE -O $SURF_DIR --enable_thirdbody --enable_falloff -E 0.005
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

quick_performance_test() {
    for m in "Hydrogen" "GRI" "A2" "C5" "JetA" "Butane" "TwoButonane" "IsoButene" "NHeptane" "IsoOctane" "Toluene" "ThreeMethylHeptane" "GasolineSurrogate" "C8_C18_Blends" "NHexadecane" "MethylFiveDeconate" "MethylDeconateNHeptane" "TwoMethylnonadecane"
    do
        echo $m
        adaptive-testing $m plug_flow_reactor -P -E 0.000001 -S PlatinumLarge -ASF -L
        # adaptive-testing $m plug_flow_reactor -L -E 0.001 -S PlatinumLarge
    done
}

run_all_surface_studies() {
    surface_threshold_study
    surface_reaction_study
    surface_performance_study
}
