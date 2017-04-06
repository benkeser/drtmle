#!/bin/bash

# Generic shell script for submitting many thousand short running jobs to a
# slurm cluster. The script goes through 4 phases of a typical HPC pipeline:
# 1. Data preparation, 2. parallel execution, 3. error correction, 4. merge
# of data. It assumes that you want to launch another script (e.g. using R
# or Python) which is set in the SCRIPT variable below. In this case it
# assumes that SCRIPT takes several command line args each of which tell
# SCRIPT to execute one of the 4 different phases. Arguments for example.R:
# Phase 0: ./example.R listsize  (get array size - not run via slurm)
# Phase 1: ./example.R prepare   (initial data preparation)
# Phase 2: ./example.R run xx    (run listsize # of jobs)
# Phase 3: ./example.R run xx    (only the failed jobs having no outfile) 
# Phase 4: ./example.R merge     (merge the output files) 
# This script (submit slurm) can take 2 arguments, SCRIPT and ANALYSIS

##################### Change these constants ##############################

ANALYSIS='newSix'      # change for every analysis you run (2nd arg)
MAILDOM='@fhcrc.org'   # your email domain (for receiving error messages)
MAXARRAYSIZE=1000          # set to 0 if you are not using slurm job arrays
MYSCRATCH="./scratch"  # location of your persistent scratch dir
PARTITION='restart'        # the queue on your cluster that allows short jobs
RESULTDIR="./out"  # This is a folder in permanent storage
SCRIPT='./centSix.R'       # your code as (R or Python) script (1st arg)
STEPSIZE=1          # number of consecutive loops in SCRIPT to run in
                           # the same job / node (increase for short jobs)
                           # MAXARRAYSIZE MUST be divisible by STEPSIZE

############## typically you don't have to change anything below here #######

username=$(id -nu)
listsize=$(${SCRIPT} listsize)   # call SCRIPT to get the number of loops

[[ -z $1 ]] || SCRIPT=$1   # get SCRIPT and ANALYSIS from command line
oldanalysis=$ANALYSIS
[[ -z $2 ]] || ANALYSIS=$2
MYSCRATCH=${MYSCRATCH//"$oldanalysis"/"$ANALYSIS"}
RESULTDIR=${RESULTDIR//"$oldanalysis"/"$ANALYSIS"}
hash ${SCRIPT} 2>/dev/null || { echo "${SCRIPT} not found, exiting !!" ; exit ; }
[[ -d ${MYSCRATCH}/out ]] && echo "${MYSCRATCH}/out already exists !!"
mkdir -p "${MYSCRATCH}/out"
mkdir -p "${MYSCRATCH}/run"
mkdir -p "${RESULTDIR}"
export MYSCRATCH           # set static vars globally so cluster jobs launched
export RESULTDIR           # from this script can read them as an environment
export STEPSIZE            # variables

echo "  using $SCRIPT with args $@..."

# run the part below if $RETRYFAILED environment var is not set or ""
if [[ -z $RETRYFAILED ]]; then 

    # Preparing input data serially, pass 'prepare' parameter to SCRIPT,
    echo "  submitting ${SCRIPT} prepare in ${MYSCRATCH}..." 
    sbatch --dependency=singleton --job-name=${ANALYSIS} \
           --partition=${PARTITION} --time=0-3 --requeue \
           --mail-type=FAIL --mail-user="${username}{MAILDOM}" \
           --output="${MYSCRATCH}/out/${ANALYSIS}.prepare.%J" \
           --wrap="${SCRIPT} prepare 0 $3 $4 $5 $6 $7 $8 $9"

    # get a colon separated list of job ids (normally a single id)
    wj=$(squeue -o "%A" -h -u ${username} -n ${ANALYSIS} -S i | tr "\n" ":")
    waitforjobs=${wj%?} #remove the last character (:) from string

    # Running multiple job arrays, pass 'run' and 'id' to each SCRIPT (1h each)
    # dividing job list into multiple job arrays
    numjobs=$((${listsize}/${MAXARRAYSIZE}+1))
    arrid=0
    arrupper=$MAXARRAYSIZE
    for i in $(seq 1 $numjobs); do
        if [[ $i -eq $numjobs ]]; then
            # set last arraysize to remainder 
            arrupper=$((${listsize}-${arrid}))
        fi
        echo "  submitting ${SCRIPT} run in ${MYSCRATCH} with "
        echo "    STEPSIZE ${STEPSIZE}, arrid ${arrid} depends on job ${waitforjobs}..."
        sbatch --dependency=afterok:${waitforjobs} --job-name=${ANALYSIS} \
               --array=1-${arrupper}:${STEPSIZE} --partition=${PARTITION} \
               --mail-type=FAIL --mail-user="${username}${MAILDOM}" --time=0-1 \
               --output="${MYSCRATCH}/out/${ANALYSIS}.run.${arrid}_%a_%A.%J" \
               --wrap="${SCRIPT} run ${arrid} $3 $4 $5 $6 $7 $8 $9" --requeue
        arrid=$((${arrid}+${MAXARRAYSIZE}))
    done

    thisscript=$0
    if ! [[ "$0" = /* ]]; then
        thisscript=$(pwd)/$0
    fi
    # submit control job that waits for parallel jobs to finish and then re-submits failed jobs (4h)
    echo "  submitting control job with args '$@' "
    echo "    that waits and resubmits failed jobs..."
    sbatch --dependency=singleton --job-name=${ANALYSIS} \
           --partition=${PARTITION} --requeue --time=0-4 \
           --mail-type=FAIL --mail-user="${username}${MAILDOM}" \
           --output="${MYSCRATCH}/out/${ANALYSIS}.correct.%J" \
           --wrap="RETRYFAILED=$numjobs ${thisscript} $1 $2 $3 $4 $5 $6 $7 $8 $9"

else
# run the part below if $RETRYFAILED environment var IS SET by previous control job
# the control job is needed because we cannot submit retry jobs right at the 
# beginning as we do not yet know which jobs will fail.

    # re-run all failed jobs where the appropriate output file is missing  (2h)
    # --requeue works only for preempted jobs 
    for i in $(seq 1 $RETRYFAILED); do
        id=$((${i}*${STEPSIZE}-${STEPSIZE}+1))
        if ! [[ -f "${MYSCRATCH}/run/${i}-run.dat" ]]; then
            echo "  re-submitting ${SCRIPT} run ${id}"
            sbatch --dependency=singleton --job-name=${ANALYSIS} \
                   --partition=${PARTITION} --requeue --time=0-2 \
                   --mail-type=FAIL --mail-user="${username}${MAILDOM}" \
                   --output="${MYSCRATCH}/out/${ANALYSIS}.run2.${i}.%J" \
                   --wrap="${SCRIPT} run ${id} $3 $4 $5 $6 $7 $8 $9"
        fi
    done

    # Merge results serially, pass 2 parameters to SCRIPT, 2nd is optional (5h)
    echo "  submitting ${SCRIPT} merge ${listsize}..."
    sbatch --dependency=singleton --job-name=${ANALYSIS} \
           --partition=${PARTITION} --requeue --time=0-5 \
           --mail-type=END,FAIL --mail-user="${username}${MAILDOM}" \
           --output="${MYSCRATCH}/out/${ANALYSIS}.merge.%J" \
           --wrap="${SCRIPT} merge ${listsize} $3 $4 $5 $6 $7 $8 $9"

fi

echo -e "  monitor output with this command:"
echo -e "  tail -f ${MYSCRATCH}/out/${ANALYSIS}.*"