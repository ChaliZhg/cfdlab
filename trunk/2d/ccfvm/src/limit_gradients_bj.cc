#include "fv.h"
#include "limiter.h"

using namespace std;

//------------------------------------------------------------------------------
// Modify gradient using Barth-Jespersen limiter
//------------------------------------------------------------------------------
void FiniteVolume::limit_gradients_bj ()
{
   vector<ConVar> Umin (grid.n_vertex);
   vector<ConVar> Umax (grid.n_vertex);

   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      Umin[i] =  1.0e20;
      Umax[i] = -1.0e20;
   }

   // For each vertex, find min and max of surrounding values
   for(unsigned int i=0; i<grid.n_cell; ++i)
      for(unsigned int j=0; j<3; ++j)
      {
         unsigned int v = grid.cell[i].vertex[j];
         
         Umin[v].min(conserved[i]);
         Umax[v].max(conserved[i]);
      }

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      phi[i] = 1.0;
      for(unsigned int j=0; j<3; ++j)
      {
         unsigned int v = grid.cell[i].vertex[j];
         
         Vector dr = grid.vertex[v].coord - grid.cell[i].centroid;
         
         double rho  = conserved[i].density    + dr * grad_rho[i];
         double rhoU = conserved[i].momentum.x + dr * grad_rhoU[i];
         double rhoV = conserved[i].momentum.y + dr * grad_rhoV[i];
         double E    = conserved[i].energy     + dr * grad_E[i];
         minmax (rho, rhoU, rhoV, E, Umin[v], Umax[v], conserved[i], phi[i]);
      }
   }

}
