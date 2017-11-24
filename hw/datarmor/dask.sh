#!/bin/csh
# Configure these values to change the size of your dask cluster
#PBS -N dask
#PBS -q mpi_2
#PBS -l select=2:ncpus=28:mpiprocs=7:ompthreads=7:mem=60g
#PBS -l walltime=24:00:00

# for infos about the queue: qstat -Qf mpi_2

# get the path right
setenv PATH ${HOME}/.miniconda3/envs/pangeo/bin:${PATH}
#setenv OMP_NUM_THREADS 1

env

# go into directory where job was launched
cd $PBS_O_WORKDIR

rm -f scheduler.json
mpirun --np 14 dask-mpi --nthreads 4 --memory-limit 1e9 --interface ib0 
# this requests 1Go per MPI process
# datarmor has 128Go per node, 128/7 = 18Go per mpi process


