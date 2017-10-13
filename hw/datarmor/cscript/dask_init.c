#include <mpi.h>
#include <stdio.h>

// Compilation:
// mpicc -o dask_init dask_init.c 

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    //printf("Hello world from processor %s, rank %d"
    //       " out of %d processors\n",
    //       processor_name, world_rank, world_size);
    //printf("%d", world_rank);

    // Print off commands to setup dask cluster
    char command[100];
    int nthreads=1;
    char *com0 ="export PATH=\"$HOME/.miniconda2/envs/daskdist/bin:$PATH\";export OMP_NUM_THREADS=1;";
    if ( world_rank == 0 ) {
      sprintf(command,"%s dask-scheduler --interface ib0 --scheduler-file $HOME/scheduler.json > %s.scheduler.out 2> %s.scheduler.err", com0, argv[1], argv[1]);
    }
    else {
      sprintf(command,"%s dask-worker --scheduler-file $HOME/scheduler.json --nthreads %d --memory-limit 0.03 --local-directory $TMPDIR --interface ib0 --name worker-%d > %s.worker.%03d.out 2> %s.worker.%03d.err;", com0, nthreads, world_rank, argv[1], world_rank, argv[1], world_rank);
    }
    system(command);
    sleep(8640000);

    // Finalize the MPI environment.
    MPI_Finalize();

}
