#!/bin/csh
# Configure these values to change the size of your dask cluster
#PBS -N dworker
#PBS -q mpi_1
#PBS -l select=1:ncpus=28:mpiprocs=7:ompthreads=7:60gb
#PBS -l walltime=24:00:00

# for infos about the queue: qstat -Qf mpi_1
# to get info about a job:   tracejob job_no    qstat -f job_no

# get the path right
setenv PATH ${HOME}/.miniconda3/envs/pangeo/bin:${PATH}
#setenv OMP_NUM_THREADS 1

# go into directory where job was launched
cd $PBS_O_WORKDIR

mpirun --np 7 dask-mpi --nthreads 4 --memory-limit 1e9 --interface ib0 --no-scheduler --local-directory ${PBS_O_WORKDIR}

# this requests 1Go per MPI process
# datarmor has 128Go per node, 128/7 = 18Go per mpi process

