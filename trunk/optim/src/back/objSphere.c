double objFun(int n, double *x)
{
   int i;
   double cost = 0.0;

   for(i=0; i<n; i++) cost += x[i]*x[i];

   return cost;

}
