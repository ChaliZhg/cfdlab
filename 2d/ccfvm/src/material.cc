#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "material.h"
#include "constants.h"

using namespace std;

extern Dimension dim;

//------------------------------------------------------------------------------
// Do some initializations
//------------------------------------------------------------------------------
void Material::initialize ()
{
   Cp = gamma * gas_const / (gamma - 1.0);

   if(model == euler)
      assert (mu_ref == 0.0);
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void Material::num_flux (const ConVar& con_left,
                         const ConVar& con_right,
                         const Vector& normal,
                         const FluxData& data,
                         Flux& flux) const
{
   PrimVar left  = con2prim(con_left);
   PrimVar right = con2prim(con_right);
   
   switch (flux_scheme)
   {
      case lxf:
         lxf_flux (left, right, normal, flux);
         break;

      case roe:
         roe_flux (left, right, normal, flux);
         break;

      case hllc:
         hllc_flux (left, right, normal, flux);
         break;
         
      case kepes_roe:
         kepes_roe_flux (left, right, normal, flux);
         break;

      case kfvs:
         kfvs_flux (left, right, normal, flux);
         break;

      default:
         cout << "num_flux: unknown flux " << flux_scheme << endl;
         exit (0);
   }
}

//------------------------------------------------------------------------------
// Flux on slip walls
//------------------------------------------------------------------------------
void Material::slip_flux (const PrimVar& state,
                          const Vector&  normal,
                          Flux&          flux) const
{
   flux.mass_flux     = 0.0;
   flux.momentum_flux = normal * state.pressure;
   flux.energy_flux   = 0.0;
}
