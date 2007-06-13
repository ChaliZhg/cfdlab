#include<stdio.h>
#include "size.h"
#include "flo.h"

void getBaseShape(int *n, int *nl, int *idl, double *xbl, double *ybl,
                  int *nu, int *idu, double *xbu, double *ybu)
{
   FILE *fp;
   int i, idummy;
   extern int myproc;

   fp = fopen("base.dat", "r");
   /* Lower surface */
   fscanf(fp, "%d", nl);
   printf("Number of points on lower surface = %d\n", *nl);
   for(i = 0; i < *nl; i++)
      fscanf(fp, "%d%lf%lf%d", &idl[i], &xbl[i], &ybl[i], &idummy);
   /* Upper surface */
   fscanf(fp, "%d", nu);
   printf("Number of points on upper surface = %d\n", *nu);
   for(i = 0; i < *nu; i++)
      fscanf(fp, "%d%lf%lf%d", &idu[i], &xbu[i], &ybu[i], &idummy);
   fclose(fp);

   *n = *nl + *nu;

   if(myproc == 0) {
      printf("Writing initial shape to shape.0\n");
      fp = fopen("shape.0", "w");
      for(i = 0; i < (*nl); i++)
         fprintf(fp, "%20.10e %20.10e %6d\n", xbl[i], ybl[i], idl[i]);
      for(i = 0; i < (*nu); i++)
         fprintf(fp, "%20.10e %20.10e %6d\n", xbu[i], ybu[i], idu[i]);
      fclose(fp);
   }

}
