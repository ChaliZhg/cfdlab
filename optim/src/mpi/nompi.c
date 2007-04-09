#include <stdio.h>
#include<mpi.h>
#include "MPI.h"

void mpi_init()
{
   myproc = 0;
   nproc  = 1;
}

/* Assign processor for each evaluation */
void mpi_assign(int n)
{
   int i;

   for(i = 1; i <= n; i++) {
      proc[i] = 0;
   }
}

/* Distribute cost function to all processess */
void mpi_distribute(int n, double *f)
{
}

/* Terminate mpi process */
void mpi_finish()
{
}
