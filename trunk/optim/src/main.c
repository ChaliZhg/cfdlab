#include<stdio.h>
#include<stdlib.h>
#include "size.h"
#include "shape.h"

int main()
{
   void getBaseShape(int *, int *, double *, double *);
   void testShapeDef();
   void multiDS_driver();

   thickness = 0.06;
   npu = npl = 5;

   getBaseShape(&nsp, idx, xb, yb);
   printf("Number of boundary points = %d\n", nsp);

   /*testShapeDef();
   exit(0);*/

   multiDS_driver();

   return 1;
}
