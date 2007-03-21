#include<math.h>

double HicksHenne(int n, double *a, double x)
{
   int i;
   double xh, tmp1, tmp2, tmp3, tmp4, y, hpower = 3.0;

   y = 0.0;
   for(i = 0; i < n; i++) {
      xh = ((double) (i + 3)) / (n + 5);
      tmp1 = log(0.5) / log(xh);
      tmp2 = M_PI * pow(x, tmp1);
      tmp3 = sin(tmp2);
      tmp4 = pow(tmp3, hpower);
      y += a[i] * tmp4;
   }

   return y;

}
