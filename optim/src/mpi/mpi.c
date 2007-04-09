#include <stdio.h>
#include <mpi.h>
#include "MPI.h"

void mpi_init()
{
   int argc, rank, size;
   char **argv;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   myproc = rank;
   nproc  = size;
   printf("Number of proc    = %d\n", nproc);
   printf("Rank of this proc = %d\n", myproc);
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
      MPI_Bcast(&f[i], 1, MPI_DOUBLE_PRECISION, proc[i], MPI_COMM_WORLD);
      printf("Vertex = %d, sender = %d, receiver = %d, cost = %f\n", i,
             proc[i], myproc, f[i]);
   }
}

/* Terminate mpi process */
void mpi_finish()
{
   MPI_Finalize();
}
