/* Deform base shape yb by applying Hicks-Henne functions */
void newShape(int npu, double *au, int npl, double *al, int nl, double *xbl,
              double *ybl, int nu, double *xbu, double *ybu, double thick,
              double *yl, double *yu)
{
   int i;
   double HicksHenne(int, double *, double);

   for(i = 0; i < nl; i++)
      yl[i] = ybl[i] + thick * HicksHenne(npl, al, xbl[i]);

   for(i = 0; i < nu; i++)
      yu[i] = ybu[i] + thick * HicksHenne(npu, au, xbu[i]);

}
