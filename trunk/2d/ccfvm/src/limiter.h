#ifndef __LIMITER_H__
#define __LIMITER_H__

#include <cmath>

//------------------------------------------------------------------------------
// minmod function
//------------------------------------------------------------------------------
inline
double minmod (const double &a, const double &b)
{
   double result;

   if( a < 0.0 && b < 0.0 )
      result = -std::min(-a, -b);
   else if( a > 0.0 && b > 0.0 )
      result = std::min(a, b);
   else
      result = 0.0;

   return result;
}

//------------------------------------------------------------------------------
// Update phi using Barth-Jespersen scheme
//------------------------------------------------------------------------------
inline
void minmax (const double& rho,
             const double& rhoU,
             const double& rhoV,
             const double& E,
             const ConVar& Umin,
             const ConVar& Umax,
             const ConVar& con,
             ConVar&       phi)
{
   // density
   if(rho > Umax.density)
   {
      double fact = (Umax.density - con.density) / 
                    (rho - con.density);
      phi.density = std::min(phi.density, fact);
   }
   else if(rho < Umin.density)
   {
      double fact = (Umin.density - con.density) / 
                    (rho - con.density);
      phi.density = std::min(phi.density, fact);
   }

   // x momentum
   if(rhoU > Umax.momentum.x)
   {
      double fact = (Umax.momentum.x - con.momentum.x) / 
                    (rhoU - con.momentum.x);
      phi.momentum.x = std::min(phi.momentum.x, fact);
   }
   else if(rhoU < Umin.momentum.x)
   {
      double fact = (Umin.momentum.x - con.momentum.x) / 
                    (rhoU - con.momentum.x);
      phi.momentum.x = std::min(phi.momentum.x, fact);
   }

   // y momentum
   if(rhoV > Umax.momentum.y)
   {
      double fact = (Umax.momentum.y - con.momentum.y) / 
                    (rhoV - con.momentum.y);
      phi.momentum.y = std::min(phi.momentum.y, fact);
   }
   else if(rhoV < Umin.momentum.y)
   {
      double fact = (Umin.momentum.y - con.momentum.y) / 
                    (rhoV - con.momentum.y);
      phi.momentum.y = std::min(phi.momentum.y, fact);
   }

   // energy
   if(E > Umax.energy)
   {
      double fact = (Umax.energy - con.energy) / 
                    (E - con.energy);
      phi.energy = std::min(phi.energy, fact);
   }
   else if(E < Umin.energy)
   {
      double fact = (Umin.energy - con.energy) / 
                    (E - con.energy);
      phi.energy = std::min(phi.energy, fact);
   }

}

#endif
