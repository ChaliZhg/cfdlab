#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

   // --------------------------------------------------------------------------
   // HLLC flux
   // Code borrowed from SU2 v2.0.2
   // --------------------------------------------------------------------------
   void Material::hllc_flux
   (
    const PrimVar  &W_l,
    const PrimVar  &W_r,
    const Vector   &normal,
    Flux           &flux
    ) const
   {
      double area = normal.norm();
      Vector unormal = normal / area;

      double rho_l_sqrt = std::sqrt(W_l.density);
      double rho_r_sqrt = std::sqrt(W_r.density);
      double fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt);
      double fact_r = 1.0 - fact_l;
      
      double v_l_unormal = W_l.velocity * unormal;
      double v_r_unormal = W_r.velocity * unormal;
      Vector velocity = W_l.velocity * fact_l + W_r.velocity * fact_r;
      double vel_unormal = velocity * unormal;
      double v2 = velocity.square();
      
      // sound speed
      double c_l = std::sqrt(gamma * W_l.pressure / W_l.density);
      double c_r = std::sqrt(gamma * W_r.pressure / W_r.density);
      
      double h_l = c_l * c_l / (gamma-1.0) + 0.5 * W_l.velocity.square();
      double h_r = c_r * c_r / (gamma-1.0) + 0.5 * W_r.velocity.square();
      
      // energy per unit mass
      double e_l = W_l.pressure/W_l.density/(gamma-1.0) + 0.5 * W_l.velocity.square();
      double e_r = W_r.pressure/W_r.density/(gamma-1.0) + 0.5 * W_r.velocity.square();
      
      // roe average
      double h = h_l * fact_l + h_r * fact_r;
      double c = std::sqrt( (gamma-1.0) * (h - 0.5*v2) );
      
      // speed of sound at l and r
      double s_l = std::min(vel_unormal-c, v_l_unormal-c_l);
      double s_r = std::max(vel_unormal+c, v_r_unormal+c_r);

      // speed of contact
      double s_m = (W_l.pressure - W_r.pressure
                    - W_l.density * v_l_unormal * (s_l-v_l_unormal)
                    + W_r.density * v_r_unormal * (s_r-v_r_unormal))
      /(W_r.density*(s_r-v_r_unormal) - W_l.density*(s_l-v_l_unormal));
      
      // Pressure at right and left (Pressure_j=Pressure_i) side of contact surface
      double pStar = W_r.density * (v_r_unormal-s_r)*(v_r_unormal-s_m) + W_r.pressure;

      if (s_m >= 0.0) {
         if (s_l > 0.0)
         {
            flux.mass_flux = W_l.density * v_l_unormal;
            flux.momentum_flux = W_l.velocity * W_l.density * v_l_unormal + unormal * W_l.pressure;
            flux.energy_flux = W_l.density * h_l * v_l_unormal;
         }
         else
         {
            double invSLmSs = 1.0/(s_l-s_m);
            double sLmuL = s_l-v_l_unormal;
            double rhoSL = W_l.density*sLmuL*invSLmSs;
            Vector rhouSL = (W_l.velocity*W_l.density*sLmuL+unormal*(pStar-W_l.pressure))*invSLmSs;
            double eSL = (sLmuL*e_l*W_l.density-W_l.pressure*v_l_unormal+pStar*s_m)*invSLmSs;
            
            flux.mass_flux = rhoSL*s_m;
            flux.momentum_flux = rhouSL*s_m + unormal*pStar;
            flux.energy_flux = (eSL+pStar)*s_m;
         }
      }
      else
      {
         if (s_r >= 0.0)
         {
            double invSRmSs = 1.0/(s_r-s_m);
            double sRmuR = s_r-v_r_unormal;
            double rhoSR = W_r.density*sRmuR*invSRmSs;
            Vector rhouSR = (W_r.velocity*W_r.density*sRmuR+unormal*(pStar-W_r.pressure))*invSRmSs;
            double eSR = (sRmuR*e_r*W_r.density-W_r.pressure*v_r_unormal+pStar*s_m)*invSRmSs;
            
            flux.mass_flux = rhoSR*s_m;
            flux.momentum_flux = rhouSR*s_m + unormal*pStar;
            flux.energy_flux = (eSR+pStar)*s_m;
         }
         else
         {
            flux.mass_flux = W_r.density*v_r_unormal;
            flux.momentum_flux = W_r.velocity*W_r.density*v_r_unormal + unormal * W_r.pressure;
            flux.energy_flux = W_r.density * h_r * v_r_unormal;
         }
      }
      
      flux *= area;
   }

