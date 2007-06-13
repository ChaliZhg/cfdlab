#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "size.h"
#include "shape.h"

void testShapeDef()
{
   int nparam;
   int i;
   FILE *fpt;
   double yl[NSPMAX], yu[NSPMAX], au[NPARMAX], al[NPARMAX];
   double dx, dy;
   void deformGrid();

   nparam = npu + npl;

   for(i = 0; i < nl; i++)
      al[i] = 0.1 * sin(i);
   for(i = 0; i < nu; i++)
      au[i] = 0.1 * cos(i);
   newShape(npu, au, npl, al, nl, xbl, ybl, nu, xbu, ybu, thickness, yl, yu);

   fpt = fopen("shape.dat", "w");
   for(i = 0; i < nl; i++)
      fprintf(fpt, "%f %f %f\n", xbl[i], ybl[i], yl[i]);
   for(i = 0; i < nu; i++)
      fprintf(fpt, "%f %f %f\n", xbu[i], ybu[i], yu[i]);
   fclose(fpt);

   /* Write shape deformation */
   fpt = fopen("P000/def.dat", "w");
   fprintf(fpt, "%d\n", nsp);
   for(i = 0; i < nl; i++) {
      dx = 0.0;
      dy = yl[i] - ybl[i];
      fprintf(fpt, "%6d %20.10e %20.10e\n", idl[i], dx, dy);
   }
   for(i = 0; i < nu; i++) {
      dx = 0.0;
      dy = yu[i] - ybu[i];
      fprintf(fpt, "%6d %20.10e %20.10e\n", idu[i], dx, dy);
   }
   fclose(fpt);

   deformGrid();
}
