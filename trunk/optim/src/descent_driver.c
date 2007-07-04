#include<stdio.h>
#include<stdlib.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

void descent_driver(int *iopt, double *ropt)
{
   int nparam;
   int i, maxiter;
   double a[NPARMAX];
   double cc, ce, lmin, lstart;
   double objFun(int, double *);
   void descent(int, double *, double, double, double, int);

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
   descent(nparam, a, cc, ce, lstart, maxiter);

}

/* Simple Steepest descent algorithm */
void descent(int nparam, double *a, double cc, double ce,
             double step, int maxiter)
{
   int i, iter = 0;
   double xold[1000], xnew[1000], gold[1000], gnew[1000], cost, tol, maxtol;

   for(i = 0; i < nparam; i++)
      xnew[i] = a[i];

   /* Main iteration loop of steepest descent */
   while(iter < maxiter && tol < maxtol) {

      for(i = 0; i < nparam; i++)
         xold[i] = xnew[i];

      /* calculate flow and cost */
      cost = objFun(nparam, xnew);

      /* calculate gradient */
      objGrad(nparam, xnew, gnew);

      /* update design variable */
      for(i = 0; i < nparam; i++)
         xnew[i] = xold[i] - step * gnew[i];

      /* check for convergence */
   }

}
