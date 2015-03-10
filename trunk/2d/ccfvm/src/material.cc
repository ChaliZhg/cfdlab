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
void Material::num_flux (const PrimVar& left0,
                         const PrimVar& right0,
                         const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         const FluxData& data,
                         Flux& flux) const
{
   switch (flux_scheme)
   {
      case lxf:
         lxf_flux (left, right, normal, flux);
         break;

      case roe:
         roe_flux (left, right, normal, flux);
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

//------------------------------------------------------------------------------
// Euler Flux Calculation
//------------------------------------------------------------------------------
void Material::euler_flux (const PrimVar& prim, 
                           const Vector&  normal,
                           Flux&          flux) const
{
   // Enthalpy 
   double h  = gamma * gas_const * prim.temperature / (gamma-1.0) + 
               0.5 * prim.velocity.square();
   // Normal velocity
   double vn = prim.velocity * normal;
   double density = Density (prim);

   flux.mass_flux = density * vn;
   flux.momentum_flux = normal * prim.pressure + 
                        prim.velocity * flux.mass_flux;
   flux.energy_flux = h * flux.mass_flux;
}

//------------------------------------------------------------------------------
// viscous flux: TODO check formulae
// This is for boundary faces only.
//------------------------------------------------------------------------------
void Material::viscous_flux (const double&  radius,
                             const bool&    adiabatic,
                             const PrimVar& state, 
                             const PrimVar& state_avg,
                             const Vector&  dU,
                             const Vector&  dV,
                             const Vector&  dW,
                             const Vector&  dT,
                             const Vector&  normal, 
                             Flux&          flux
                             ) const
{
   // state_avg is used for material coefficients
   double mu = viscosity (state_avg.temperature);
   double k = mu * Cp / prandtl;

   // Heat flux
   double q;
   if(adiabatic)
      q = 0.0;
   else
      q = -k * (dT * normal);

   // Divergence of velocity
   double div = dU.x + dV.y;

   // Extra term in case of axisymmetric
   // x is radial direction
   double stx = 0, sty = 0;
   if(dim == axi)
   {
      div += state.velocity.x / radius;

      stx = mu * (dW.x - state.velocity.z / radius);
      sty = mu * dW.y;
   }

   // Shear stress tensor: symmetric, compute only upper part
   double sxx = 2.0 * mu * (dU.x - (1.0/3.0) * div);
   double syy = 2.0 * mu * (dV.y - (1.0/3.0) * div);
   double sxy = mu * (dU.y + dV.x);

   flux.mass_flux = 0.0;
   flux.momentum_flux.x = -(sxx * normal.x + sxy * normal.y);
   flux.momentum_flux.y = -(sxy * normal.x + syy * normal.y);
   flux.momentum_flux.z = -(stx * normal.x + sty * normal.y);
   flux.energy_flux = flux.momentum_flux * state.velocity + q;

}

//------------------------------------------------------------------------------
// viscous flux: TODO check formulae
// This is for interior faces.
//------------------------------------------------------------------------------
void Material::viscous_flux (const double&  radius,
                             const PrimVar& state, 
                             const Vector&  dU,
                             const Vector&  dV,
                             const Vector&  dW,
                             const Vector&  dT,
                             const Vector&  normal0, Flux& flux0,
                             const Vector&  normal1, Flux& flux1,
                             const Vector&  normal2, Flux& flux2
                             ) const
{
   double mu = viscosity (state.temperature);
   double k = mu * Cp / prandtl;

   // Heat flux
   double q0 = -k * (dT * normal0);
   double q1 = -k * (dT * normal1);
   double q2 = -k * (dT * normal2);

   // Divergence of velocity
   double div = dU.x + dV.y;

   // Extra term in case of axisymmetric
   // x is radial direction
   double stx = 0, sty = 0;
   if(dim == axi)
   {
      div += state.velocity.x / radius;

      stx = mu * (dW.x - state.velocity.z / radius);
      sty = mu * dW.y;
   }

   // Shear stress tensor: symmetric, compute only upper part
   double sxx = 2.0 * mu * (dU.x - (1.0/3.0) * div);
   double syy = 2.0 * mu * (dV.y - (1.0/3.0) * div);
   double sxy = mu * (dU.y + dV.x);

   flux0.mass_flux       = 0.0;
   flux0.momentum_flux.x = -(sxx * normal0.x + sxy * normal0.y);
   flux0.momentum_flux.y = -(sxy * normal0.x + syy * normal0.y);
   flux0.momentum_flux.z = -(stx * normal0.x + sty * normal0.y);
   flux0.energy_flux     = flux0.momentum_flux * state.velocity + q0;

   flux1.mass_flux       = 0.0;
   flux1.momentum_flux.x = -(sxx * normal1.x + sxy * normal1.y);
   flux1.momentum_flux.y = -(sxy * normal1.x + syy * normal1.y);
   flux1.momentum_flux.z = -(stx * normal1.x + sty * normal1.y);
   flux1.energy_flux     = flux1.momentum_flux * state.velocity + q1;

   flux2.mass_flux       = 0.0;
   flux2.momentum_flux.x = -(sxx * normal2.x + sxy * normal2.y);
   flux2.momentum_flux.y = -(sxy * normal2.x + syy * normal2.y);
   flux2.momentum_flux.z = -(stx * normal2.x + sty * normal2.y);
   flux2.energy_flux     = flux2.momentum_flux * state.velocity + q2;
}

//------------------------------------------------------------------------------
// Source terms in axisymmetric flow
//------------------------------------------------------------------------------
void Material::axisymmetric_source(const double&  radius,
                                   const PrimVar& state,
                                   const Vector&  dU, 
                                   const Vector&  dV, 
                                   const Vector&  dW, 
                                   Flux&          source) const
{
   source.mass_flux   = 0.0;
   source.energy_flux = 0.0;

   double mu = viscosity (state.temperature);
   double density = Density(state);
   double div = dU.x + dV.y + state.velocity.x / radius;
   double stt = 2.0 * mu * ( state.velocity.x / radius - (1.0/3.0) * div );
   double stx = mu * (dW.x - state.velocity.z / radius);

   // radial equation
   source.momentum_flux.x = - state.pressure - density * pow(state.velocity.z,2.0) + stt;

   // axial equation
   source.momentum_flux.y = 0.0;

   // tangential equation
   source.momentum_flux.z = density * state.velocity.x * state.velocity.z - stx;
}
