#include <vector>
#include "fv.h"

using namespace std;

// Jameson's residual smoothing
// Uses Jacobi iterations
void FiniteVolume::smooth_residual()
{
   if(!param.smooth_res) return;

   static const double eps = 0.5;
   static const unsigned int niter = 2;

   for(unsigned int iter=0; iter<niter; ++iter)
   {
      if(iter==0)
         for(unsigned int i=0; i<grid.n_vertex; ++i)
            residual1[i] = residual[i];
      else
         for(unsigned int i=0; i<grid.n_vertex; ++i)
            residual1[i] = residual2[i];

      for(unsigned int i=0; i<grid.n_vertex; ++i)
      {
         residual2[i] = residual[i];
         unsigned int nnbr = grid.vertex[i].nbr_vertex.size();
         for(unsigned int j=0; j<nnbr; ++j)
         {
            unsigned int nbr = grid.vertex[i].nbr_vertex[j];
            residual2[i] += residual1[nbr] * eps;
         }
         residual2[i] *= 1.0/(1.0 + nnbr*eps);
      }
   }

   // copy smoothed residual
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      residual[i] = residual2[i];

}
