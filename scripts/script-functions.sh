#!/bin/bash
# create directorys
check_args_and_dirs() {
    if [ -z "$1" ]
    then
        echo "No args file given"
        exit 1
    else
        export ADD_ARGS="$(<$1)"
        export DIRNAME=${ADD_ARGS##* }
        if [ ! -d $DIRNAME ]
        then
            mkdir $DIRNAME
        fi
    fi
    # make slurm-output dir
    if [ ! -d "slurm-output" ]
    then
        mkdir slurm-output
    fi
}
# check and set threshold bounds
check_set_threshold_bounds() {
    # starting threshold
    if [ -z "$TSTART" ]
    then
        export TSTART=0
    fi
    # ending threshold
    if [ -z "$TEND" ]
    then
        export TEND=18
    fi
}

# function to run approximate preconditioned single jobs
approximate_precon_single() {
    # set threshold bounds
    check_set_threshold_bounds
    # loop through thresholds
    for th in $(seq $TSTART $TEND)
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL -L -v -M -P -T $THR --prefix approx $ADD_ARGS"
        echo $JOB_OPTIONS
        echo
        # run the job with set options
        sbatch -J "$CURR_MODEL-approx-single" ./batches/jobs-single.sh --mem=$AMS
    done
}

# function to run approximate preconditioned mpi jobs
approximate_precon_mpi() {
    # set threshold bounds
    check_set_threshold_bounds
    # loop through thresholds
    for th in $(seq $TSTART $TEND)
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL -L -v -M -P -T $THR --prefix approx $ADD_ARGS"
        echo $JOB_OPTIONS
        echo
        # run the job with set options
        sbatch -J "$CURR_MODEL-approx-mpi" ./batches/jobs-mpi.sh --mem=$AMS
    done
}

# function to run fully analytical preconditioned single jobs
analytical_precon_single() {
    # set threshold bounds
    check_set_threshold_bounds
    # loop through thresholds
    for th in $(seq $TSTART $TEND)
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL -L -v -M -P -T $THR --prefix analyt --skip_thirdbody --skip_falloff --analyt_temp_derivs $ADD_ARGS"
        echo $JOB_OPTIONS
        echo
        # run the job with set options
        sbatch -J "$CURR_MODEL-analyt-single" ./batches/jobs-single.sh --mem=$AMS
    done
}

# function to run fully analytical preconditioned mpi jobs
analytical_precon_mpi() {
    # set threshold bounds
    check_set_threshold_bounds
    # loop through thresholds
    for th in $(seq $TSTART $TEND)
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL -L -v -M -P -T $THR --prefix analyt --skip_thirdbody --skip_falloff --analyt_temp_derivs $ADD_ARGS"
        echo $JOB_OPTIONS
        echo
        # run the job with set options
        sbatch -J "$CURR_MODEL-analyt-mpi" ./batches/jobs-mpi.sh --mem=$AMS
    done
}

# function to run fully analytical preconditioned single jobs
mass_single() {
    export JOB_OPTIONS="$CURR_MODEL -L -v $ADD_ARGS"
    echo $JOB_OPTIONS
    echo
    # run the job with set options
    sbatch -J "$CURR_MODEL-mass-mpi" ./batches/jobs-single.sh --mem=$AMS
}

# function to run fully analytical preconditioned mpi jobs
mass_mpi() {
    export JOB_OPTIONS="$CURR_MODEL -L -v $ADD_ARGS"
    echo $JOB_OPTIONS
    echo
    # run the job with set options
    sbatch -J "$CURR_MODEL-mass-mpi" ./batches/jobs-mpi.sh --mem=$AMS
}

# function to run fully analytical preconditioned single jobs
moles_single() {
    export JOB_OPTIONS="$CURR_MODEL -L -v -M $ADD_ARGS"
    echo $JOB_OPTIONS
    echo
    # run the job with set options
    sbatch -J "$CURR_MODEL-moles-mpi" ./batches/jobs-single.sh --mem=$AMS
}

# function to run fully analytical preconditioned mpi jobs
moles_mpi() {
    export JOB_OPTIONS="$CURR_MODEL -L -v -M $ADD_ARGS"
    echo $JOB_OPTIONS
    echo
    # run the job with set options
    sbatch -J "$CURR_MODEL-moles-mpi" ./batches/jobs-mpi.sh --mem=$AMS
}

skip_mass() {
    export SKIP_MASS=1
}

skip_moles() {
    export SKIP_MOLES=1
}

skip_approx() {
    export SKIP_APPROX=1
}

skip_analyt() {
    export SKIP_ANALYT=1
}

# Define single jobs
define_runners() {
    if [ -z "$1" ]
    then
        echo "define runners not given runner type"
        return
    fi
    # Approximate preconditioner
    if [ -z "$SKIP_APPROX" ]
    then
        RUNNERS+=("approximate_precon_$1")
    fi
    # Analytical preconditioner
    if [ -z "$SKIP_ANALYT" ]
    then
        RUNNERS+=("analytical_precon_$1")
    fi
    # Mass run
    if [ -z "$SKIP_MASS" ]
    then
        RUNNERS+=("mass_$1")
    fi
    # Moles run
    if [ -z "$SKIP_MOLES" ]
    then
        RUNNERS+=("moles_$1")
    fi
}

# TODO: current not functional
profile_preconditioning() {
    source ~/.functions
    export CURR_DIR=$(pwd)
    export EX_DIR='../examples'
    export CANTERA_DIR='../../cantera/'
    # Rebuild cantera
    cd $CANTERA_DIR
    scons build -j 12
    scons install prefix=$CONDA_PREFIX
    # Build example
    cd $CURR_DIR
    cd $EX_DIR
    scons;
    ./dev-test.out;
    # profile example
    gprof -Q -b dev-test.out gmon.out>$CURR_DIR/profiles/$1
    # -Q no call graph
    # -D ignore non-functions
}
