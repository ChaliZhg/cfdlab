/* Torczon's multi-directionsl search
 * n = number of design parameters
 * x = double array of design parameters; initial value on input
 * lmin = minimum edge length of simplex for convergence criterion
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void multiDS(int n, double *x, double cc, double ce, double lmin,
             double lstart, int maxiter)
{
   int i, imin, replaced, iter = 0;
   double **xs, **xr, **xe, **xc, *fs, *fr, *fe, *fc, fsmin, frmin,
      femin, fcmin, ssize;
   void initSimplex(int, double *, double **, double);
   void printSimplex(int, int, double **, double *);
   void findBest(int, double **, double *, int *, double *);
   double objFun(int, double *);
   void copySimplex(int, double **, double **, double *, double *);
   double simplexSize(int, double **);
   void vecAdd(int, double *, double *, double *, double);
   double dmin(int, double *);

   /* Initial size of simplex */
   ssize = lstart;

   /* Check validity of input parameters */
   if(cc <= 0.0 || cc >= 1.0) {
      printf("multiDS: contraction coefficient must be in (0,1)\n");
      exit(0);
   }

   if(ce <= 1.0) {
      printf("multiDS: expandion coefficient must be > 1\n");
      exit(0);
   }

   if(ssize < lmin) {
      printf("multiDS: starting simplex size is < minimum\n");
      printf("         give lstart > lmin\n");
      exit(0);
   }

   printf("Parameters for search:\n");
   printf("   Contraction factor     = %e\n", cc);
   printf("   Expansion   factor     = %e\n", ce);
   printf("   Starting simplex size  = %e\n", ssize);
   printf("   Minimum  simplex size  = %e\n", lmin);
   printf("   Maximum number of iter = %d\n", maxiter);

   /* Allocate memory */
   xs = (double **) calloc((n + 1), sizeof(double *));
   xr = (double **) calloc((n + 1), sizeof(double *));
   xe = (double **) calloc((n + 1), sizeof(double *));
   xc = (double **) calloc((n + 1), sizeof(double *));
   fs = (double *) calloc(n + 1, sizeof(double));
   fr = (double *) calloc(n + 1, sizeof(double));
   fe = (double *) calloc(n + 1, sizeof(double));
   fc = (double *) calloc(n + 1, sizeof(double));
   for(i = 0; i < n + 1; i++) {
      xs[i] = (double *) calloc(n, sizeof(double));
      xr[i] = (double *) calloc(n, sizeof(double));
      xe[i] = (double *) calloc(n, sizeof(double));
      xc[i] = (double *) calloc(n, sizeof(double));
   }


   /* Initialize the simplex */
   initSimplex(n, x, xs, ssize);

   /* Calculate initial function values */
   printf("Initial simplex and function values:\n");
   for(i = 0; i < n + 1; i++) {
      fs[i] = objFun(n, xs[i]);
   }
   printSimplex(0, n, xs, fs);

   /* Find best vertex and put in first position */
   findBest(n, xs, fs, &imin, &fsmin);

   /* Main iteration loop */
   while(ssize > lmin && iter < maxiter) {
      printf("Iteration = %d\n\n", iter + 1);

      replaced = 0;
      while(!replaced && ssize > lmin) { /* inner repeat loop */

         /* rotation step */
         printf("   Rotation:\n");
         for(i = 1; i <= n; i++) {
            vecAdd(n, xs[0], xs[i], xr[i], 1.0);
            fr[i] = objFun(n, xr[i]);
         }
         printSimplex(1, n, xr, fr);

         frmin = dmin(n, fr);
         replaced = (frmin < fs[0]) ? 1 : 0;
         if(replaced) {
            /* expansion step */
            printf("   Expansion:\n");
            for(i = 1; i <= n; i++) {
               vecAdd(n, xs[0], xs[i], xe[i], ce);
               fe[i] = objFun(n, xe[i]);
            }
            printSimplex(1, n, xe, fe);

            femin = dmin(n, fe);
            if(femin < frmin)
               copySimplex(n, xe, xs, fe, fs); //accept expansion
            else
               copySimplex(n, xr, xs, fr, fs); //accept rotation
         }
         else {
            /* contraction step */
            printf("   Contraction step:\n");
            for(i = 1; i <= n; i++) {
               vecAdd(n, xs[0], xs[i], xc[i], -cc);
               fc[i] = objFun(n, xc[i]);
            }
            printSimplex(1, n, xc, fc);

            fcmin = dmin(n, fc);
            replaced = (fcmin < fs[0]) ? 1 : 0;
            copySimplex(n, xc, xs, fc, fs); //accept contraction
         }

         /* Length of smallest edge in simplex */
         ssize = simplexSize(n, xs);

      }                         /* End of inner repeat loop */

      ++iter;

      /* Find best vertex and put in first position */
      findBest(n, xs, fs, &imin, &fsmin);

      printf("\n");
      printf("Minimum length of simplex = %12.4e\n", ssize);
      printf("Minimum function value    = %12.4e\n", fs[0]);
      printf("-------------------------------------------------\n");
   }                            /* End of main iteration loop */

   /* Copy best vertex for output */
   for(i = 0; i < n; i++)
      x[i] = xs[0][i];

   /* Best vertex found */
   printf("Best vertex:\n");
   for(i = 0; i < n; i++) printf("%e ", x[i]);
   printf("\n");


   /* Free memory */
   for(i = 0; i < n + 1; i++) {
      free(xs[i]);
      free(xr[i]);
      free(xe[i]);
      free(xc[i]);
   }
   free(xs);
   free(xr);
   free(xe);
   free(xc);
   free(fs);
   free(fr);
   free(fe);
   free(fc);
}

/* Initialize the simplex */
void initSimplex(int n, double *x, double **xs, double r)
{
   int i, j;

   /* Set every vertex to initial vertex */
   for(i = 0; i < n + 1; i++)
      for(j = 0; j < n; j++)
         xs[i][j] = x[j];

   for(i = 1; i < n + 1; i++)
      xs[i][i - 1] = xs[i][i - 1] + r;

}

/* Find length of smallest edge in simplex */
double simplexSize(int n, double **x)
{
   int i, j, k;
   double dx, rmin, r;

   rmin = 1.0e20;
   for(i = 0; i <= n; i++)
      for(j = i + 1; j <= n; j++) {
         r = 0.0;
         for(k = 0; k < n; k++) {
            dx = x[i][k] - x[j][k];
            r += dx * dx;
         }
         r = sqrt(r);
         rmin = (r < rmin) ? r : rmin;
      }

   return rmin;
}

/* Find best vertex and put it in first position */
void findBest(int n, double **xs, double *fs, int *imin, double *fsmin)
{
   int i, j, IMIN;
   double FSMIN, ftmp, *xt;

   xt = (double *) calloc(n, sizeof(double));

   /* Find best vertex */
   IMIN = 0;
   FSMIN = fs[0];
   for(i = 1; i < n + 1; i++)
      if(fs[i] < FSMIN) {
         FSMIN = fs[i];
         IMIN = i;
      }

   /* Put best vertex in 0'th position */
   for(j = 0; j < n; j++)
      xt[j] = xs[0][j];
   for(j = 0; j < n; j++)
      xs[0][j] = xs[IMIN][j];
   for(j = 0; j < n; j++)
      xs[IMIN][j] = xt[j];
   ftmp = fs[0];
   fs[0] = fs[IMIN];
   fs[IMIN] = ftmp;

   *imin = IMIN;
   *fsmin = FSMIN;

   free(xt);
}

/* Print simplex and function values */
void printSimplex(int i1, int n, double **x, double *f)
{
   int i, j;
   for(i = i1; i <= n; i++) {
      printf("     ");
      for(j = 0; j < n; j++)
         printf("%12.4e", x[i][j]);
      printf("%12.4e\n", f[i]);
   }
}

/* x = x0 - c*(xi - x0) */
void vecAdd(int n, double *x0, double *xi, double *x, double c)
{
   int j;

   for(j = 0; j < n; j++)
      x[j] = x0[j] - c * (xi[j] - x0[j]);
}

/* Find minimum of f[1], f[2], ..., f[n] */
double dmin(int n, double *f)
{
   int i;
   double fm;

   fm = f[1];
   for(i = 2; i <= n; i++)
      fm = (f[i] < fm) ? f[i] : fm;
   return fm;
}

/* x2[i][j] = x1[i][j], i=1,...,n, j=0,...,n-1 */
void copySimplex(int n, double **x1, double **x2, double *f1, double *f2)
{
   int i, j;

   for(i = 1; i <= n; i++) {
      f2[i] = f1[i];
      for(j = 0; j < n; j++)
         x2[i][j] = x1[i][j];
   }
}
