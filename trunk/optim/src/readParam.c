#include <stdio.h>

void readParam(int *iopt, double *ropt)
{
   FILE *fpt;
   char dummy[100];

   fpt = fopen("optim.inp", "r");
   fscanf(fpt, "%s%d%d", dummy, &iopt[0], &iopt[1]);
   fscanf(fpt, "%s%d", dummy, &iopt[2]);

   fscanf(fpt, "%s%lf", dummy, &ropt[0]);
   fscanf(fpt, "%s%lf", dummy, &ropt[1]);
   fscanf(fpt, "%s%lf", dummy, &ropt[2]);
   fscanf(fpt, "%s%lf", dummy, &ropt[3]);
   fscanf(fpt, "%s%lf", dummy, &ropt[4]);
   fclose(fpt);

}
