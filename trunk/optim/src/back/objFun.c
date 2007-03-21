#include<stdio.h>
#include "size.h"
#include "shape.h"

double objFun(int n, double *x)
{
   int i, iter;
   double au[NPARMAX], al[NPARMAX], y[NSPMAX];
   double cost = 0.0, dx, dy;
   double residue, cl, cd;
   FILE *fpt;
   void deformGrid(void);
   void runSolver(void);

   /* Separate Hicks-Henne params for lower and upper surface */
   for(i = 0; i < npl; i++)
      al[i] = x[i];
   for(i = 0; i < npu; i++)
      au[i] = x[i + npl];

   /* Generate new shape */
   newShape(npu, au, npl, al, nsp, xb, yb, thickness, y);

   /* Write shape deformation */
   fpt = fopen("def.dat", "w");
   fprintf(fpt, "%d\n", nsp);
   for(i = 0; i < nsp; i++) {
      dx = 0.0;
      dy = y[i] - yb[i];
      fprintf(fpt, "%6d %20.10e %20.10e\n", idx[i] + 1, dx, dy);
   }
   fclose(fpt);

   /* Deform the grid */
   deformGrid();

   /* Run flow solver */
   runSolver();

   /* Read flo solution */
   fpt = fopen("FLO.OUT", "r");
   fscanf(fpt,"%d%lf%lf%lf", &iter, &residue, &cl, &cd);
   fclose(ftp);

   return cost;

}
