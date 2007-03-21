#include<stdio.h>
#include<stdlib.h>
#include "size.h"
#include "shape.h"

void testShapeDef()
{
   int nparam;
   int i;
   FILE *fpt;
   double y[NSPMAX], au[NPARMAX], al[NPARMAX];
   double dx, dy;
   void deformGrid();

   nparam = npu + npl;

   au[0] = 0.1;
   au[1] = 0.2;
   au[2] = -0.1;
   au[3] = 0.2;
   au[4] = 0.1;
   al[0] = 0.1;
   al[1] = 0.2;
   al[2] = -0.1;
   al[3] = 0.2;
   al[4] = 0.1;
   newShape(npu, au, npl, al, nsp, xb, yb, thickness, y);

   fpt = fopen("shape.dat", "w");
   for(i = 0; i < nsp; i++)
      fprintf(fpt, "%f %f %f\n", xb[i], yb[i], y[i]);
   fclose(fpt);

   /* Write shape deformation */
   fpt = fopen("def.dat", "w");
   fprintf(fpt, "%d\n", nsp);
   for(i = 0; i < nsp; i++) {
      dx = 0.0;
      dy = y[i] - yb[i];
      fprintf(fpt, "%6d %20.10e %20.10e\n", idx[i] + 1, dx, dy);
   }
   fclose(fpt);

   deformGrid();
}
