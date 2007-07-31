#include "double.h"

#include <stdio.h>
#include <stdlib.h>
#include "routines.h"
#include "param.h"

void residu(int nc, int nf, REALX * a, REALQ ** q, REALA ** res)
{
   int i, j, nl, nr;

   for(i = 0; i < nc; i++)
      for(j = 0; j < 3; j++)
         res[i][j] = 0.0;

   flux_in(a[0], res[0]);
   for(i = 1; i < nf - 1; i++) {
      nl = i - 1;
      nr = i;
      if(flux1 == 1)
         flux_ausm(a[i], q[nl], q[nr], res[nl], res[nr]);
//    else if(flux1 == 2)
//       flux_kfvs(a[i], q[nl], q[nr], res[nl], res[nr]);
      else if(flux1 == 3)
         flux_lf(a[i], q[nl], q[nr], res[nl], res[nr]);
      else {
         printf("Flux type is not known, %d\n", flux1);
         exit(0);
      }
   }
   flux_out(a[nf - 1], q[nf - 2], q[nf - 3], res[nf - 2]);

   for(i = 0; i < nc; i++)
      source_term(a[i], a[i + 1], q[i], res[i]);

}
