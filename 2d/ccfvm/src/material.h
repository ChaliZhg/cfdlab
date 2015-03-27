#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include "vec.h"
#include "primvar.h"
#include "flux.h"
#include "convar.h"

//------------------------------------------------------------------------------
// Logarithmic average: (a - b)/log(a/b)
// Numerically stable alogorithm taken from Ismail and Roe
//------------------------------------------------------------------------------
inline
double logavg(double a, double b)
{
   double xi = b/a;
   double f = (xi - 1.0) / (xi + 1.0);
   double u = f * f;

   double F;
   if (u < 1.0e-2)
   {
      double u2 = u * u;
      double u3 = u2 * u;
      F = 1.0 + u/3.0 + u2/5.0 + u3/7.0;
   }
   else
      F = log(xi)/2.0/f;

   return 0.5*(a+b)/F;
}

//------------------------------------------------------------------------------
struct FluxData
{
   double ssw;
   double ducros;
};

//------------------------------------------------------------------------------
// Material class
//------------------------------------------------------------------------------
class Material
{
   public:
      Material ()
      {
         mu_ref = 0.0;
      };
      double gamma;
      double gas_const;
      double prandtl;
      double Cp;
      double T_0, T_ref, mu_ref; // constants for sutherland law
      double omega; // exponent in power law viscosity
      enum FlowModel {euler, ns};
      FlowModel model;
      enum FluxScheme { kep, lxf, roe, hllc, kfvs, kepes, kepes_roe, kepes_roe2, kepes_rus,
                        kepes_hyb };
      FluxScheme flux_scheme;

      enum MuModel {mu_constant, mu_sutherland, mu_power};
      MuModel mu_model;

      void initialize ();
      ConVar  prim2con (const PrimVar& prim_var) const;
      PrimVar con2prim (const ConVar&  con_var) const;
      void num_flux(const ConVar&   left,
                    const ConVar&   right,
                    const Vector&   normal,
                    const FluxData& data,
                    Flux& flux) const;
      void    lxf_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void    roe_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void    roe2_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void   hllc_flux (const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
   void   kepes_roe_flux (const PrimVar& left,
                     const PrimVar& right,
                     const Vector& normal,
                     Flux& flux) const;
      void kfvs_split_flux (const double   sign,
                            const Vector&  normal,
                            const PrimVar& state,
                            Flux&          flux) const;
      void    kfvs_flux(const PrimVar& left, 
                        const PrimVar& right, 
                        const Vector& normal, 
                        Flux& flux) const;
      void    slip_flux (const PrimVar& state, 
                         const Vector& normal, 
                         Flux& flux) const;
      void    euler_flux (const PrimVar& prim, 
                          const Vector&  normal,
                          Flux& flux) const;
      void viscous_flux (const double&  radius,
                         const bool&    adiabatic,
                         const PrimVar& state, 
                         const PrimVar& state_avg,
                         const Vector&  dU, 
                         const Vector&  dV, 
                         const Vector&  dW, 
                         const Vector&  dT, 
                         const Vector&  normal, Flux& flux
                         ) const;
      void viscous_flux (const double&  radius,
                         const PrimVar& state, 
                         const Vector&  dU, 
                         const Vector&  dV, 
                         const Vector&  dW, 
                         const Vector&  dT, 
                         const Vector&  normal0, Flux& flux0,
                         const Vector&  normal1, Flux& flux1,
                         const Vector&  normal2, Flux& flux2
                         ) const;
      void axisymmetric_source(const double&  radius,
                               const PrimVar& state,
                               const Vector&  dU, 
                               const Vector&  dV, 
                               const Vector&  dW, 
                               Flux&          source) const;
      double viscosity (const double T) const;
      double total_energy (const PrimVar& state) const;
      double sound_speed (const PrimVar& state) const;
      double Mach (const PrimVar& state) const;
      double Pressure(const ConVar& state) const;

};

//------------------------------------------------------------------------------
// Convert primitive to conserved
//------------------------------------------------------------------------------
inline
ConVar Material::prim2con(const PrimVar& prim_var) const
{
   ConVar con_var;

   con_var.density  = prim_var.density;
   con_var.momentum = prim_var.velocity * prim_var.density;
   con_var.energy   = prim_var.pressure/(gamma - 1.0) +
                        0.5 * prim_var.velocity.square() * prim_var.density;

   return con_var;
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
inline
PrimVar Material::con2prim (const ConVar& con_var) const
{
   PrimVar prim_var;

   prim_var.velocity = con_var.momentum / con_var.density;
   prim_var.pressure = (gamma - 1.0) * 
        ( con_var.energy - 0.5 * con_var.momentum.square() / con_var.density );
   prim_var.density = con_var.density;

   return prim_var;
}

//------------------------------------------------------------------------------
// Viscosity coefficient according to sutherland law
//------------------------------------------------------------------------------
inline
double Material::viscosity (const double T) const
{
   switch (mu_model)
   {
      case mu_constant:
         return mu_ref;

      case mu_sutherland:
         return mu_ref * std::pow(T/T_ref, 1.5) * (T_ref + T_0) / (T + T_0);

      case mu_power:
         return mu_ref * std::pow(T/T_ref, omega);

      default:
         std::cout << "viscosity: unknown model " << mu_model << std::endl;
         abort ();
   }
}

//------------------------------------------------------------------------------
// Total energy per unit volume
//------------------------------------------------------------------------------
inline
double Material::total_energy (const PrimVar& state) const
{
   return state.pressure / (gamma - 1.0) +
          0.5 * state.density * state.velocity.square();
}

//------------------------------------------------------------------------------
// sound speed
//------------------------------------------------------------------------------
inline
double Material::sound_speed (const PrimVar& state) const
{
   return sqrt(gamma * state.pressure / state.density);
}

//------------------------------------------------------------------------------
// Mach number
//------------------------------------------------------------------------------
inline
double Material::Mach (const PrimVar& state) const
{
   double sonic = sound_speed (state);
   return state.velocity.norm() / sonic;
}

//------------------------------------------------------------------------------
// Pressure
//------------------------------------------------------------------------------
inline
double Material::Pressure(const ConVar& state) const
{
   return (gamma - 1.0) * (state.energy - 0.5 * state.momentum.square() / state.density);
}


#endif
