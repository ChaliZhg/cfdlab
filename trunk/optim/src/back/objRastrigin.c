#include<math.h>

double objFun(int n, double *x)
{
   int i;
   double cost = 10.0*n;

   for(i=0; i<n; i++) cost += x[i]*x[i] - 10.0*cos(2.0*M_PI*x[i]);

   return cost;

}
