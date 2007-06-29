#include <stdio.h>
#include <mpi.h>
#include "MPI.h"
#include "opt.h"

void mpi_init(int argc, char **argv)
{
   int rank, size, len;
   char command[100];

   printf("Initializing mpi ...\n");

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   myproc = rank;
   nproc  = size;
   printf("Number of proc    = %d\n", nproc);
   printf("Rank of this proc = %d\n", myproc);

   if(rank < 10)               sprintf(rundir, "P00%d", rank);
   if(rank>=10 && rank <100)   sprintf(rundir, "P0%d", rank);
   if(rank>=100 && rank <1000) sprintf(rundir, "P%d", rank);

   sprintf(deform, "./run.sh deform %s", rundir);
   sprintf(flosol, "./run.sh solve %s", rundir);
   sprintf(adjsol, "./run.sh adjoint %s", rundir);
   sprintf(adjmesh, "./run.sh adjmesh %s", rundir);

   printf("Run directory  = %s\n", rundir);
   printf("Deform command = %s\n", deform);
   printf("Flow   command = %s\n", flosol);

   sprintf(command,"rm -rf %s && mkdir %s", rundir, rundir);
   system(command);
}

/* Assign processor for each evaluation */
void mpi_assign(int n)
{
   int i, count;

   count = 0;
   proc[0] = -1;
   for(i = 1; i <= n; i++) {
      proc[i] = count;
      count = count + 1;
      if(count == nproc)
         count = 0;
      printf("Evaluation = %d, Processor = %d\n", i, proc[i]);
   }
}

/* Distribute cost function to all processess */
void mpi_distribute(int n, double *f)
{
   int i;

   for(i = 1; i <= n; i++) {
      MPI_Bcast(&f[i], 1, MPI_DOUBLE, proc[i], MPI_COMM_WORLD);
      printf("Vertex = %d, sender = %d, receiver = %d, cost = %f\n", i,
             proc[i], myproc, f[i]);
   }
}

/* Terminate mpi process */
void mpi_finish()
{
   MPI_Finalize();
}
