#include<stdio.h>
#include "size.h"
#include "flo.h"

void getBaseShape(int *n, int *idx, double *x, double *yb)
{
   FILE *fp;
   int i, np, nt, idummy, ptype[NPMAX];
   double xg[NPMAX], yg[NPMAX];

   fp = fopen("mesh.0", "r");
   fscanf(fp, "%d%d", &np, &nt);
   printf("Number of points    in grid = %d\n", np);
   printf("Number of triangles in grid = %d\n", nt);
   for(i = 0; i < np; i++)
      fscanf(fp, "%d%lf%lf%d", &idummy, &xg[i], &yg[i], &ptype[i]);
   fclose(fp);

   *n = 0;
   for(i = 0; i < np; i++)
      if(ptype[i] == SOLID) {
         x[*n] = xg[i];
         yb[*n] = yg[i];
         idx[*n] = i;
         ++(*n);
      }

   printf("Writing initial shape to shape.0\n");
   fp = fopen("shape.0", "w");
   for(i = 0; i < (*n); i++)
      fprintf(fp, "%20.10e %20.10e %6d\n", x[i], yb[i], idx[i]);
   fclose(fp);

}
