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
      
      // Apply positivity limiter. Make pressure at face centers to be positive.
      const double eps = 1.0e-13;
      for(unsigned int j=0; j<3; ++j)
      {
         unsigned int f = grid.cell[i].face[j];
         Vector dr = grid.face[f].centroid - grid.cell[i].centroid;
         
         ConVar con;
         con.density    = conserved[i].density    + (dr * grad_rho[i] ) * phi[i].density;
         con.momentum.x = conserved[i].momentum.x + (dr * grad_rhoU[i]) * phi[i].momentum.x;
         con.momentum.y = conserved[i].momentum.y + (dr * grad_rhoV[i]) * phi[i].momentum.y;
         con.energy     = conserved[i].energy     + (dr * grad_E[i]   ) * phi[i].energy;
         
         double pre = param.material.Pressure (con);
         if(pre < eps)
         {
            //std::cout << "Pressure is negative\n";
            //exit(0);
            phi[i] = 0.0;
         }
      }
   }

}
