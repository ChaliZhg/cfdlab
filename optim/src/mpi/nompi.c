#include <stdio.h>
#include <mpi.h>
#include "MPI.h"
#include "opt.h"

void mpi_init(int argc, char **argv)
{
   //extern char rundir[], deform[], flosol[], adjsol[], adjmesh[];
   char command[100];

   myproc = 0;
   nproc  = 1;

   sprintf(rundir, "P000");
   sprintf(deform, "./run.sh deform P000");
   sprintf(flosol, "./run.sh solve P000");
   sprintf(adjsol, "./run.sh adjoint P000");
   sprintf(adjmesh, "./run.sh adjmesh P000");

   printf("Run directory = %s\n", rundir);
   printf("Deform command = %s\n", deform);
   printf("Flow   command = %s\n", flosol);

   sprintf(command,"rm -rf %s && mkdir %s", rundir, rundir);
   system(command);
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
