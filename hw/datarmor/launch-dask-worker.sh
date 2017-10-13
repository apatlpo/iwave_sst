#!/bin/csh
# Configure these values to change the size of your dask cluster
#PBS -q mpi_8
#PBS -l walltime=24:00:00

# for infos about the queue: qstat -Qf mpi_8

# get the path for mpirun
#source /usr/share/Modules/3.2.10/init/csh
# module purge
 #module load mpt
# module load impi/5.1.3.258
 #module load intel-cmkl-17/17.0.2.174
 #module load   NETCDF/4.3.3.1-mpt-intel2016  # faster (9:03.32)  but result not exact  as caparmor
 #module load   NETCDF/4.3.3.1-mpt-intel2015  # result ok slower(10:49.27 )
 module list
#setenv mpiproc `cat $PBS_NODEFILE  |wc -l`
setenv PATH ${HOME}/.miniconda2/envs/sst/bin:${PATH}
setenv OMP_NUM_THREADS 1

# go into directory where job was launched
cd $PBS_O_WORKDIR

echo $MPI_LAUNCH  # not there anymore
cat $PBS_NODEFILE | uniq >  "$PBS_JOBID.workers.nodefile"
setenv SCHEDULER ${HOME}/scheduler.json

#time $MPI_LAUNCH -n 8  dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
#time $MPI_LAUNCH -np 8  dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
time mpirun -np 8  dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
    --no-scheduler --local-directory $TMPDIR \
    --scheduler-file=$SCHEDULER


# --np 6 dask-mpi --nthreads 6 --memory-limit 22e9 --interface ib0 \
#    --no-scheduler --local-directory /glade/scratch/$USER \
#    --scheduler-file=$SCHEDULER


#!/bin/bash
#PBS -N dask-worker
#PBS -q economy
#PBS -A UCLB0022
#PBS -l select=1:ncpus=36:mpiprocs=6:ompthreads=6
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -m abe

# Setup Environment
#module purge
#source activate pangeo

# Setup dask worker
#SCHEDULER=$HOME/scheduler.json
#mpirun --np 6 dask-mpi --nthreads 14 --memory-limit 22e9 --interface ib0 \
#    --no-scheduler --local-directory $TMPDIR \
#    --scheduler-file=$SCHEDULER
