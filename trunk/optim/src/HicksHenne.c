#include<math.h>

double HicksHenne(int n, double *a, double x)
{
   int i;
   double y;
   double HicksHenneFun(int, int, double);

   y = 0.0;
   for(i = 0; i < n; i++)
      y += a[i] * HicksHenneFun(i, n, x);

   return y;

}

/* i'th Hicks-Henne function */
double HicksHenneFun(int i, int n, double x)
{
   double xh, tmp1, tmp2, tmp3, tmp4, hpower = 3.0;

   xh = ((double) (i + 3)) / (n + 5);
   tmp1 = log(0.5) / log(xh);
   tmp2 = M_PI * pow(x, tmp1);
   tmp3 = sin(tmp2);
   tmp4 = pow(tmp3, hpower);
   return tmp4;
}
