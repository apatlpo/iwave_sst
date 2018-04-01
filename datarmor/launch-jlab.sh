#!/bin/bash
set -e

# Usage:
#   bash
#   source activate pangeo
#   ./launch-jlab.sh 

source activate pangeon

# create a directory to store temporary dask data
#mkdir -p $DATAWORK/dask 
#rm -rf $DATAWORK/dask/*
#
#mkdir -p $SCRATCH/dask 
#echo "Clean up ${SCRATCH}/dask"
#rm -rf $SCRATCH/dask/*

echo "Launching job"
s=`qsub jlab.pbs`
sjob=${s%.*}
echo ${s}

#qstat -f ${sjob}

# block until the scheduler job starts
while true; do
    status=`qstat ${sjob} | tail -n 1`
    #echo ${status}
    if [[ ${status} =~ " R " ]]; then
        break
    fi
    sleep 1
done

#./jlab.py ${s}

#

#default=$HOME
#notebook=${2:-$default}
#echo "Setting up Jupyter Lab, Notebook dir: ${notebook}"
#./setup-jlab.py --log_level=DEBUG --jlab_port=8877 --dash_port=8878 \
#    --notebook_dir $notebook --scheduler_file $SCHEDULER


