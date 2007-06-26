#include <stdio.h>
#include "opt.h"

void readParam(int *iopt, double *ropt)
{
   FILE *fpt;
   char dummy[100];
   int nsolid, i;
   double rdummy;

   fpt = fopen("optim.inp", "r");
   fscanf(fpt, "%s%d%d", dummy, &iopt[0], &iopt[1]);
   fscanf(fpt, "%s%d", dummy, &iopt[2]);
   fscanf(fpt, "%s%d", dummy, &iopt[3]);

   fscanf(fpt, "%s%lf", dummy, &ropt[0]);
   fscanf(fpt, "%s%lf", dummy, &ropt[1]);
   fscanf(fpt, "%s%lf", dummy, &ropt[2]);
   fscanf(fpt, "%s%lf", dummy, &ropt[3]);
   fscanf(fpt, "%s%lf", dummy, &ropt[4]);
   fclose(fpt);

   /* If costfun==1, i.e., pressure matching problem, then read the target
    * pressure distribution. Note that the order of points in this file must
    * match the order in which the flow-solver outputs the pressure
    * coefficient into WALL.DAT */
   if(iopt[3] == 1) {
      printf("Reading cp0.dat ... ");
      fpt = fopen("cp0.dat", "r");
      fscanf(fpt, "%d", &nsolid);
      printf("there are %5d points in there\n", nsolid);
      for(i = 0; i < nsolid; i++)
         fscanf(fpt, "%lf%lf%lf", &rdummy, &CP0[i], &rdummy);
      fclose(fpt);
   }

}
