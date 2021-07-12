#!/bin/bash

date=`date +%Y-%m-%d-T%H-%M-%S`
memory=12g
jobname="calibration_$date"
num_cores=1
walltime=12:00:00
row=1

while getopts m:r:T:c:j: option
do
    case "${option}"
        in
        m) memory=${OPTARG};;
        r) row=${OPTARG};;
        T) walltime=${OPTARG};;
        c) num_cores=${OPTARG};;
        j) jobname=${OPTARG}
    esac
done

usage() {
    echo "
    usage: ./subCal.sh [-m memory] [-r row] [-T walltime] [-c num_cores] [-j jobname]
    Starts a profound calibration simulation

    options:
    -m memory         per node, as #[k|m|g] (default: $memory)
    -r row            start of batch
    -T walltime       as hh:mm:ss, max compute time (default $walltime)
    -c num_cores	  How many cores to request and run the job on (default: $num_cores)
    -j jobname	      name of analysis for organization (default: calibration_date)
    "
    exit 0
}
echo starting $jobname in $HOME/scratch/profound
mkdir -p $HOME/scratch/profound
mkdir $HOME/scratch/profound/$jobname
sbatch --output=$HOME/scratch/profound/$jobname/slurm.out -J $jobname -t $walltime --mem=$memory -c $num_cores ./calibration.R $row $HOME/scratch/profound $num_cores

