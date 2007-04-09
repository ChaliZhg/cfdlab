#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "size.h"
#include "shape.h"

int main()
{
   void getBaseShape(int *, int *, double *, double *);
   void testShapeDef();
   void multiDS_driver();

   mpi_init();

   thickness = 0.06;
   npu = npl = 1;

   /* getBaseShape(&nsp, idx, xb, yb);
   printf("Number of boundary points = %d\n", nsp);
   */

   /*testShapeDef();
      exit(0); */

   multiDS_driver();

   mpi_finish();
   return 1;
}
