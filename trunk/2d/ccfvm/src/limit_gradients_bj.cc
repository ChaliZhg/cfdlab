#include "fv.h"
#include "limiter.h"

using namespace std;

//------------------------------------------------------------------------------
// Modify gradient using Barth-Jespersen limiter
//------------------------------------------------------------------------------
void FiniteVolume::limit_gradients_bj ()
{
   vector<PrimVar> pmin (grid.n_vertex);
   vector<PrimVar> pmax (grid.n_vertex);

   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      pmin[i] = primitive[i];
      pmax[i] = primitive[i];
      phi[i]  = 1.0;
   }

   // For each vertex, find min and max of surrounding values
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      unsigned int v0 = grid.face[i].vertex[0];
      unsigned int v1 = grid.face[i].vertex[1];
      pmin[v0].min(primitive[v1]);
      pmax[v0].max(primitive[v1]);
      pmin[v1].min(primitive[v0]);
      pmax[v1].max(primitive[v0]);
   }

   // Compute limiter
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      unsigned int n0 = grid.face[i].vertex[0];
      unsigned int n1 = grid.face[i].vertex[1];
      Vector dr = grid.vertex[n1].coord - grid.vertex[n0].coord;

      double T0 = primitive[n0].density     + 0.5 * (dT[n0] * dr);
      double u0 = primitive[n0].velocity.x  + 0.5 * (dU[n0] * dr);
      double v0 = primitive[n0].velocity.y  + 0.5 * (dV[n0] * dr);
      double w0 = primitive[n0].velocity.z  + 0.5 * (dW[n0] * dr);
      double p0 = primitive[n0].pressure    + 0.5 * (dP[n0] * dr);

      minmax (T0, u0, v0, w0, p0, pmin[n0], pmax[n0], primitive[n0], phi[n0]);

      double T1 = primitive[n1].density     - 0.5 * (dT[n1] * dr);
      double u1 = primitive[n1].velocity.x  - 0.5 * (dU[n1] * dr);
      double v1 = primitive[n1].velocity.y  - 0.5 * (dV[n1] * dr);
      double w1 = primitive[n1].velocity.z  - 0.5 * (dW[n1] * dr);
      double p1 = primitive[n1].pressure    - 0.5 * (dP[n1] * dr);

      minmax (T1, u1, v1, w1, p1, pmin[n1], pmax[n1], primitive[n1], phi[n1]);
   }

}
