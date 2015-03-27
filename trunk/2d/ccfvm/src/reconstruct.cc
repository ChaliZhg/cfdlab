#include <cmath>
#include "fv.h"
#include "limiter.h"

#define EPSILON  1.0e-14
#define limit_albada(a,b)  max(0.0, (2*a*b + EPSILON)/(a*a + b*b + EPSILON))

using namespace std;

//------------------------------------------------------------------------------
// First order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_first
(
 const unsigned int& f,
 vector<ConVar>&    state
) const
{
   // Left state
   int cl = grid.face[f].lcell;
   state[0] = conserved[cl];

   int cr = grid.face[f].rcell;
   if(cr >= 0)
      state[1] = conserved[cr];
   else
      state[1] = state[0];
}

//------------------------------------------------------------------------------
// Second order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_second
(
 const unsigned int& f,
 const Vector& p,
 vector<ConVar>&    state
) const
{
   int cl = grid.face[f].lcell;
   int cr = grid.face[f].rcell;

   Vector  drl = p - grid.cell[cl].centroid;

   // left state
   state[0].density    = conserved[cl].density    + drl * grad_rho[cl];
   state[0].momentum.x = conserved[cl].momentum.x + drl * grad_rhoU[cl];
   state[0].momentum.y = conserved[cl].momentum.y + drl * grad_rhoV[cl];
   state[0].energy     = conserved[cl].energy     + drl * grad_E[cl];

   if(cr >= 0)
   {
      Vector  drr = p - grid.cell[cr].centroid;
      
      // right state
      state[1].density    = conserved[cr].density    + drr * grad_rho[cr];
      state[1].momentum.x = conserved[cr].momentum.x + drr * grad_rhoU[cr];
      state[1].momentum.y = conserved[cr].momentum.y + drr * grad_rhoV[cr];
      state[1].energy     = conserved[cr].energy     + drr * grad_E[cr];
   }
   else
      state[1] = state[0];
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_minmax
(
 const unsigned int& f,
 const Vector& p,
 vector<ConVar>&    state
 ) const
{
   int cl = grid.face[f].lcell;
   int cr = grid.face[f].rcell;
   
   Vector  drl = p - grid.cell[cl].centroid;
   
   // left state
   state[0].density    = conserved[cl].density    + (drl * grad_rho[cl] ) * phi[cl].density;
   state[0].momentum.x = conserved[cl].momentum.x + (drl * grad_rhoU[cl]) * phi[cl].momentum.x;
   state[0].momentum.y = conserved[cl].momentum.y + (drl * grad_rhoV[cl]) * phi[cl].momentum.y;
   state[0].energy     = conserved[cl].energy     + (drl * grad_E[cl]   ) * phi[cl].energy;
   
   if(cr >= 0)
   {
      Vector  drr = p - grid.cell[cr].centroid;
      
      // right state
      state[1].density    = conserved[cr].density    + (drr * grad_rho[cr] ) * phi[cr].density;
      state[1].momentum.x = conserved[cr].momentum.x + (drr * grad_rhoU[cr]) * phi[cr].momentum.x;
      state[1].momentum.y = conserved[cr].momentum.y + (drr * grad_rhoV[cr]) * phi[cr].momentum.y;
      state[1].energy     = conserved[cr].energy     + (drr * grad_E[cr]   ) * phi[cr].energy;
   }
   else
      state[1] = state[0];
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct (const unsigned int& f,
                                const Vector& p,
                                vector<ConVar>&    state) const
{
   switch(param.reconstruct_scheme)
   {
      // First order
      case Parameter::first:
         reconstruct_first (f, state);
         break;
         
      // Second order
      case Parameter::second:
         reconstruct_second (f, p, state);
         break;

      // Barth-Jespersen or MinMax scheme
      case Parameter::bj:
         reconstruct_minmax (f, p, state);
         break;

      default:
         cout << "reconstruct: unknown reconstruction scheme = " 
              << param.reconstruct_scheme << endl;
         exit (0);
   }

   // Modify velocity by Thornber scheme for low mach case
   //prec_thornber(state);
}
