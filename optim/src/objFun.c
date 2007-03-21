#include<stdio.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

double objFun(int n, double *x)
{
   int i, iter;
   double au[NPARMAX], al[NPARMAX], y[NSPMAX];
   double cost = 0.0, dx, dy;
   double residue, tmp, liftpenalty;
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
   fscanf(fpt, "%d%lf%lf%lf", &iter, &residue, &cl, &cd);
   fclose(fpt);
   printf("Number of iterations = %d\n", iter);
   printf("Flow residue         = %e\n", residue);
   printf("Lift coefficient, cl = %e\n", cl);
   printf("Drag coefficient, cd = %e\n", cd);

   tmp = 1.0 - cl / clref;
   tmp = (tmp > 0.0) ? tmp : 0.0;
   liftpenalty = 1.0e4 * tmp;
   cost = cd / cdref + liftpenalty;

   return cost;

}
