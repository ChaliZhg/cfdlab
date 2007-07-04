#include<stdio.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

/* Find gradient of cost function */
void objGrad(int n, double *x, double *grad)
{
   int i, j, p, iter, indx[NSPMAX], indx1[NSPMAX];
   double au[NPARMAX], al[NPARMAX];
   double xcoord[NSPMAX], ycoord[NSPMAX], gx[NSPMAX], gy[NSPMAX], gal[NSPMAX],
      gau[NSPMAX];
   FILE *fpt;
   char gradfile[100];
   void runAdjoint(void);
   void runMeshAdjoint(void);

   sprintf(gradfile, "%s/DEF.OUT", rundir);

   /* Separate Hicks-Henne params for lower and upper surface */
   for(i = 0; i < npl; i++)
      al[i] = x[i];
   for(i = 0; i < npu; i++)
      au[i] = x[i + npl];

   /* Run adjoint flow solver */
   runAdjoint();

   /* Run adjoint mesh solver */
   runMeshAdjoint();

   /* Read gradient wrt shape points */
   fpt = fopen(gradfile, "r");
   for(i = 0; i < nsp; i++)
      fscanf(fpt, "%d%lf%lf%lf%lf", &indx[i], &xcoord[i], &ycoord[i], &gx[i],
             &gy[i]);
   fclose(fpt);

   /* Identify points on lower and upper surface */
   for(i = 0; i < nsp; i++) {
      indx1[i] = 0;
      for(j = 0; j < nl; j++)
         if(indx[i] == idl[j])
            indx1[i] = -1;
      for(j = 0; j < nu; j++)
         if(indx[i] == idu[j])
            indx1[i] = +1;
      if(indx1[i] == 0) {
         printf("objGrad: error in locating point %d\n", i);
         printf("         cant decide if it is on lower or upper surface\n");
         exit(0);
      }
   }

   /* Grad wrt upper surface Hicks-Henne parameter */
   for(i = 0; i < npl; i++) {
      gal[i] = 0.0;
      for(p = 0; p < nsp; p++)
         if(indx1[p] == -1)
            gal[i] += thickness * gy[p] * HicksHenneFun(i, npl, xcoord[p]);
   }

   /* Grad wrt lower surface Hicks-Henne parameter */
   for(i = 0; i < npu; i++) {
      gau[i] = 0.0;
      for(p = 0; p < nsp; p++)
         if(indx1[p] == +1)
            gau[i] += thickness * gy[p] * HicksHenneFun(i, npu, xcoord[p]);
   }

   /* Put the gradients into single output array */
   for(i = 0; i < npl; i++)
      grad[i] = gal[i];
   for(i = 0; i < npu; i++)
      grad[i + npl] = gau[i];

}
