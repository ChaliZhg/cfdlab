#include<stdio.h>
#include "size.h"
#include "shape.h"
#include "opt.h"

/* Find gradient of cost function */
void objGrad(int n, double *x, double *grad)
{
   int i, j, p, iter, indx[NSPMAX];
   double au[NPARMAX], al[NPARMAX];
   double xcoord[NSPMAX], ycoord[NSPMAX], gx[NSPMAX], gy[NSPMAX], gal[NSPMAX],
      gau[NSPMAX];
   FILE *fpt;
   void runAdjoint(void);
   void runMeshAdjoint(void);

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
   fpt = fopen("DEF.OUT", "r");
   for(i = 0; i < nsp; i++)
      fscanf(fpt, "%d%lf%lf%lf%lf", &indx[i], &xcoord[i], &ycoord[i], &gx[i],
             &gy[i]);
   fclose(fpt);

   /* Identify points on lower and upper surface */
   for(i = 0; i < nsp; i++) {
      for(j = 0; j < npl; j++)
         if(indx[i] == idl[j])
            indx[i] = -1;
      for(j = 0; j < npu; j++)
         if(indx[i] == idu[j])
            indx[i] = +1;
   }

   /* Grad wrt upper surface Hicks-Henne parameter */
   for(i = 0; i < npl; i++) {
      gal[i] = 0.0;
      for(p = 0; p < nsp; p++)
         if(indx[p] == -1)
            gal[i] += thickness * gy[p] * HicksHenneFun(i, npl, xcoord[p]);
   }

   /* Grad wrt lower surface Hicks-Henne parameter */
   for(i = 0; i < npu; i++) {
      gau[i] = 0.0;
      for(p = 0; p < nsp; p++)
         if(indx[p] == +1)
            gau[i] += thickness * gy[p] * HicksHenneFun(i, npu, xcoord[p]);
   }

   /* Put the gradients into single output array */
   for(i = 0; i < npl; i++)
      grad[i] = gal[i];
   for(i = 0; i < npu; i++)
      grad[i + npl] = gau[i];

}
