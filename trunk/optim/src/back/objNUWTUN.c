#include<stdio.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

double objNUWTUN(int n, double *x)
{
   int i, iter;
   double au[NPARMAX], al[NPARMAX], yl[NSPMAX], yu[NSPMAX];
   double cost = 0.0, dx, dy, rdummy;
   double residue, tmp, liftpenalty;
   char deffile[100], flofile[100], cpfile[100];
   FILE *fpt;
   void deformGrid(void);
   void runSolver(void);

   sprintf(deffile, "%s/hicks.in", rundir);
   sprintf(flofile, "%s/FLO.OUT", rundir);
   sprintf(cpfile, "%s/WALL.DAT", rundir);

   printf("Hicks-Henne parameters:\n");
   for(i = 0; i < n; i++)
      printf("%5d %12.4e\n", i, x[i]);

   /* Rudimentary handling of bound constraint. We should put this in an
    * input file  */
   for(i = 0; i < n; i++)
      if(x[i] > 0.25 || x[i] < -0.25) {
         printf("Out of bounds; setting large value\n");
         return 1.0e20;
      }

   /* Write Hicks-Henne parameters */
   fpt = fopen(deffile, "w");
   for(i = 0; i < n; i++) {
      fprintf(fpt, "%25.15e\n", x[i]);
   }
   fclose(fpt);

   /* Deform the grid */
   deformGrid();

   /* Run flow solver */
   runSolver();

   /* Read flo solution */
   fpt = fopen(flofile, "r");
   fscanf(fpt, "%d%lf%lf%lf", &iter, &residue, &cl, &cd);
   fclose(fpt);
   printf("Number of iterations = %d\n", iter);
   printf("Flow residue         = %e\n", residue);
   printf("Lift coefficient, cl = %e\n", cl);
   printf("Drag coefficient, cd = %e\n", cd);

   /* Read pressure coefficient on the airfoil */
   fpt = fopen(cpfile, "r");
   for(i = 0; i < nsp; i++)
      fscanf(fpt, "%lf%lf%lf", &rdummy, &CP[i], &rdummy);
   fclose(fpt);

   /* Compute cost function */
   switch (costfun) {
      case 1:                  /* pressure matching */
         cost = 0.0;
         for(i = 0; i < nsp; i++)
            cost += 0.5 * (CP[i] - CP0[i]) * (CP[i] - CP0[i]);
         break;
      case 2:                  /* drag with lift penalty */
         tmp = 1.0 - cl / clref;
         tmp = (tmp > 0.0) ? tmp : 0.0;
         liftpenalty = 1.0e4 * tmp;
         cost = cd / cdref + liftpenalty;
         break;
      default:
         cost = 0.0;
         printf("Cost function type %4d is unknown !!!\n", costfun);
   }

   printf("Cost function        = %e\n", cost);

   return cost;

}
