#!/bin/bash
# check and wait so as not to exceed slurm job limit
slurm_job_wait() {
    n_jobs=$(squeue -u walkanth | wc -l)
    stime=1
    while [ $n_jobs -ge 100 ]
    do
        echo "Number of jobs over 100, waiting $stime seconds..."
        sleep $stime
        n_jobs=$(squeue -u walkanth | wc -l)
        # adjust sleep time
        if [ $stime -lt 10 ]
        then
            ((stime=stime+1))
        fi
    done
    echo "Currently queued jobs: $n_jobs"
}
# create directorys
check_args_and_dirs() {
    if [ -z "$1" ]
    then
        echo "No args file given"
        exit 1
    fi
    export PROBLEMS="$(<$1)"
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
    for (( th = ${TSTART}; th <= ${TEND}; th++ ))
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P -T $THR $ADD_ARGS"
        echo $JOB_OPTIONS $AMS
        echo
        # run the job with set options
        if [ -z "$SKIP_SBATCH" ]
        then
            slurm_job_wait
            sbatch -J "$CURR_MODEL-approx-single" --mem=$AMS ./batches/jobs-single.sh
        fi
        sleep $SLEEP_TIMER
    done
}

# function to run approximate preconditioned mpi jobs
approximate_precon_mpi() {
    # set threshold bounds
    check_set_threshold_bounds
    # loop through thresholds
    for (( th = ${TSTART}; th <= ${TEND}; th++ ))
    do
        if [ $th -ne 0 ]
        then
            THR="1e-$th"
        else
            THR="0"
        fi
        export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P -T $THR $ADD_ARGS"
        echo $JOB_OPTIONS $AMS
        echo
        # run the job with set options
        if [ -z "$SKIP_SBATCH" ]
        then
            slurm_job_wait
            sbatch -J "$CURR_MODEL-approx-mpi" --mem=$AMS ./batches/jobs-mpi.sh
        fi
        sleep $SLEEP_TIMER
    done
}

# function to run fully analytical preconditioned single jobs
analytical_precon_single() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P -T 0 --prefix analyt --skip_thirdbody --skip_falloff --analyt_temp_derivs $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-analyt-single" --mem=$AMS ./batches/jobs-single.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run fully analytical preconditioned mpi jobs
analytical_precon_mpi() {
    # loop through thresholds
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P -T 0 --prefix analyt --skip_thirdbody --skip_falloff --analyt_temp_derivs $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-analyt-mpi" --mem=$AMS  ./batches/jobs-mpi.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run fully analytical preconditioned single jobs
mass_single() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-mass-single" --mem=$AMS ./batches/jobs-single.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run fully analytical preconditioned mpi jobs
mass_mpi() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-mass-mpi" --mem=$AMS ./batches/jobs-mpi.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run fully analytical preconditioned single jobs
moles_single() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-moles-single" --mem=$AMS ./batches/jobs-single.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run fully analytical preconditioned mpi jobs
moles_mpi() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-moles-mpi" --mem=$AMS ./batches/jobs-mpi.sh
    fi
    sleep $SLEEP_TIMER
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

skip_sbatch() {
    export SKIP_SBATCH=1
}

reset_skips() {
    unset SKIP_ANALYT
    unset SKIP_APPROX
    unset SKIP_MASS
    unset SKIP_MOLES
    unset SKIP_SBATCH
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

job_wait_loop() {
    n_jobs=$(squeue -u $USER --format="%.18i %.40j" | grep -E $1 | wc -l)
    stime=1
    while [ $n_jobs -gt 0 ]
    do
        echo "Currently queued jobs: $n_jobs. Waiting $stime seconds..."
        sleep $stime
        n_jobs=$(squeue -u $USER --format="%.18i %.40j" | grep -E $1 | wc -l)
        # adjust sleep time
        if [ $stime -lt 10 ]
        then
            ((stime=stime+1))
        fi
    done

    sbatch -J "JOBS_ALL_COMPLETED" ./batches/completed.sh
}
