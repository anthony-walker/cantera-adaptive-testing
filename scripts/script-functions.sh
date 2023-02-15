#!/bin/bash

# function to print jobs in a specific format
slurm_job_print() {
    squeue -u $USER --format="%.18i %.40j"
}

# check and wait so as not to exceed slurm job limit
slurm_job_wait() {
    # njobs setter
    if [ -z "$NJOBS_LIMIT" ]
    then
        export NJOBS_LIMIT=1000
    fi
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        n_jobs=$(squeue -u walkanth | wc -l)
        stime=1
        while [ $n_jobs -ge $NJOBS_LIMIT ]
        do
            echo "Number of jobs over $NJOBS_LIMIT, waiting $stime seconds..."
            sleep $stime
            n_jobs=$(squeue -u walkanth | wc -l)
            # adjust sleep time
            if [ $stime -lt 10 ]
            then
                ((stime=stime+1))
            fi
        done
        echo "Currently queued jobs: $n_jobs"
    else
        sleep $SLEEP_TIMER
    fi
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

# function to run approximate preconditioned single jobs
flexible_precon_single() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P --flexible $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-flex-single" --mem=$AMS ./batches/jobs-single.sh
    fi
    sleep $SLEEP_TIMER
}

# function to run approximate preconditioned mpi jobs
flexible_precon_mpi() {
    export JOB_OPTIONS="$CURR_MODEL $PROBLEMS -v -M -P --flexible $ADD_ARGS"
    echo $JOB_OPTIONS $AMS
    echo
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        slurm_job_wait
        sbatch -J "$CURR_MODEL-flex-mpi" --mem=$AMS ./batches/jobs-mpi.sh
    fi
    sleep $SLEEP_TIMER
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

skip_flex() {
    export SKIP_FLEX=1
}

skip_steady_state() {
    export SKIP_STEADY_STATE=1
}


reset_skips() {
    unset SKIP_ANALYT
    unset SKIP_APPROX
    unset SKIP_MASS
    unset SKIP_MOLES
    unset SKIP_SBATCH
    unset SKIP_FLEX
    unset SKIP_STEADY_STATE
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
    # Moles run
    if [ -z "$SKIP_FLEX" ]
    then
        RUNNERS+=("flexible_precon_$1")
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
    export CURR_DIR=$(pwd)
    export EX_DIR='../examples'
    # export CANTERA_DIR='../../cantera/'
    # # Rebuild cantera
    # cd $CANTERA_DIR
    # scons build -j 12
    # scons install prefix=$CONDA_PREFIX
    # # Build example
    # cd $CURR_DIR
    cd $EX_DIR
    scons;
    ./dev-test.out;
    # profile example
    gprof -Q -b dev-test.out gmon.out>$CURR_DIR/profiles/$1
    # -Q no call graph
    # -D ignore non-functions
}

steady_state_wait() {
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
        steady_args="-R steady -MTS 1e-3 -MS 1e9 -O $SURF_DIR"
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
        if [ $stime -lt 60 ]
        then
            ((stime=stime+1))
        fi
    done
    python -c "from postprocess_surface import *; combine_surf_yamls()"
    sbatch -J "JOBS_ALL_COMPLETED" ./batches/completed.sh
}
