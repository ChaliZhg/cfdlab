#include <vector>
#include "fv.h"

using namespace std;

// Compute ducros switch
void FiniteVolume::compute_ducros()
{
   if(param.ducros == false)
   {
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         ducros[i] = 1;
      return;
   }

   // Loop over interior faces and accumulate flux
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      double div  = dU[i].x + dV[i].y;
      double vor  = dV[i].x - dU[i].y;

      div = div*div;
      vor = vor*vor;

      ducros[i] = min( div/(div + vor + 1.0e-6), 1.0);
   }
}
