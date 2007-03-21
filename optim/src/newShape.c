/* Deform base shape yb by applying Hicks-Henne functions */
void newShape(int npu, double *au, int npl, double *al, int nsp, double *x,
              double *yb, double thick, double *y)
{
   int i;
   double HicksHenne(int, double *, double);

   for(i = 0; i < nsp; i++)
      if(yb[i] >= 0.0)
         y[i] = yb[i] + thick * HicksHenne(npu, au, x[i]);
      else
         y[i] = yb[i] + thick * HicksHenne(npl, al, x[i]);

}
