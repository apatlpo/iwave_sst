#!/bin/csh
# Configure these values to change the size of your dask cluster
#PBS -N dask
#PBS -q mpi_2
#PBS -l select=2:ncpus=28:mpiprocs=7:ompthreads=7
#PBS -l walltime=24:00:00

# for infos about the queue: qstat -Qf mpi_8

# get the path right
setenv PATH ${HOME}/.miniconda3/envs/pangeo/bin:${PATH}
#setenv OMP_NUM_THREADS 1

# go into directory where job was launched
cd $PBS_O_WORKDIR

rm -f scheduler.json
mpirun --np 14 dask-mpi --nthreads 4 --memory-limit 24e9 --interface ib0 
# this requests 24Go per MPI process
# datarmor has 128Go per node, 128/7 = 18Go per mpi process

#cat $PBS_NODEFILE | uniq >  "$PBS_JOBID.scheduler.nodefile"
#setenv SCHEDULER ${HOME}/scheduler.json
#rm -f $SCHEDULER

#time mpirun -np 8  dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
#    --local-directory $TMPDIR --scheduler-file=$SCHEDULER

