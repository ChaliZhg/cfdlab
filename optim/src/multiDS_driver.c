#include<stdio.h>
#include<stdlib.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

void multiDS_driver(int *iopt, double *ropt)
{
   int nparam;
   int i, maxiter;
   double a[NPARMAX];
   double cc, ce, lmin, lstart;
   void multiDS(int, double *, double, double, double, double, int);

   system("rm -f flo.log");

   nparam = npu + npl;
   printf("Number of design parameters = %d\n", nparam);

   cc = ropt[1];
   ce = ropt[2];
   lmin = ropt[3];
   lstart = ropt[4];
   maxiter = iopt[2];
   costfun = iopt[3];

   /* Give some non-zero value */
   clref = cdref = 1.0;

   /* Solve for initial shape */
   for(i = 0; i < nparam; i++)
      a[i] = 0.0;
   cost0 = objFun(nparam, a);

   /* Set reference values */
   clref = cl;
   cdref = cd;

   printf("Reference values cl, cd = %e %e\n", cl, cd);

   /* Call multi-directional search */
   multiDS(nparam, a, cc, ce, lmin, lstart, maxiter);

}
