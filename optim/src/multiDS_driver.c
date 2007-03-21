#include<stdio.h>
#include<stdlib.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

void multiDS_driver()
{
   int nparam;
   int i, maxiter;
   double a[NPARMAX];
   double cc, ce, lmin, lstart, cost;
   double objFun(int, double *);
   void multiDS(int, double *, double, double, double, double, int);

   nparam = npu + npl;
   printf("Number of design parameters = %d\n", nparam);

   cc = 0.5;
   ce = 1.5;
   lmin = 1.0e-10;
   lstart = 0.5;
   maxiter = 100;

   /* Solve for initial shape */
   for(i = 0; i < nparam; i++)
      a[i] = 0.0;
   cost = objFun(nparam, a);

   /* Set reference values */
   clref = cl;
   cdref = cd;

   printf("Reference values cl, cd = %e %e\n", cl, cd);
   exit(0);

   /* Call multi-directional search */
   multiDS(nparam, a, cc, ce, lmin, lstart, maxiter);

}
