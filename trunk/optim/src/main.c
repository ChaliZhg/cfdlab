#include <stdio.h>
#include <stdlib.h>
#include "size.h"
#include "shape.h"
#include "opt.h"
#include "flo.h"
//#include "mpi/MPI.h"

int main(int argc, char **argv)
{
   void readParam(int *, double *);
   void getBaseShape(int *n, int *nl, int *idl, double *xbl, double *ybl,
                     int *nu, int *idu, double *xbu, double *ybu);
   void testShapeDef();
   void multiDS_driver(int *, double *);
   void mpi_init(int, char **);
   void mpi_finish();
   int iopt[100];
   double ropt[100];

   mpi_init(argc, argv);

   /*mpi_finish();
      return 1;
      exit(0); */

   /* Read input parameters */
   readParam(iopt, ropt);

   thickness = ropt[0];
   npu = iopt[0];
   npl = iopt[1];

   getBaseShape(&nsp, &nl, idl, xbl, ybl, &nu, idu, xbu, ybu);
   printf("Number of boundary points = %d\n", nsp);

   /*testShapeDef();
      exit(0); */

   multiDS_driver(iopt, ropt);

   mpi_finish();
   return 1;
}
