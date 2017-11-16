#!/bin/csh
# Configure these values to change the size of your dask cluster
#PBS -N dask-scheduler
#PBS -q mpi_8
#PBS -l walltime=24:00:00

# for infos about the queue: qstat -Qf mpi_8

Update this file by following dask.sh

# get the path right
setenv PATH ${HOME}/.miniconda2/envs/sst/bin:${PATH}
setenv OMP_NUM_THREADS 1

# go into directory where job was launched
cd $PBS_O_WORKDIR

cat $PBS_NODEFILE | uniq >  "$PBS_JOBID.scheduler.nodefile"
setenv SCHEDULER ${HOME}/scheduler.json
rm -f $SCHEDULER

time mpirun -np 8  dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
    --local-directory $TMPDIR --scheduler-file=$SCHEDULER

