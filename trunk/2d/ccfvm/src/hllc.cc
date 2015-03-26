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
   
   double rho_l_sqrt = sqrt(W_l.density);
   double rho_r_sqrt = sqrt(W_r.density);
   double fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt);
   double fact_r = 1.0 - fact_l;
   
   double v_l_normal = W_l.velocity * unormal;
   double v_r_normal = W_r.velocity * unormal;
   Vector velocity = W_l.velocity * fact_l + W_r.velocity * fact_r;
   double vel_normal = velocity * unormal;
   
   // sound speed
   double c_l = sqrt(gamma * W_l.pressure / W_l.density);
   double c_r = sqrt(gamma * W_r.pressure / W_r.density);
   
   double h_l = c_l * c_l / (gamma-1.0) + 0.5 * W_l.velocity.square();
   double h_r = c_r * c_r / (gamma-1.0) + 0.5 * W_r.velocity.square();
   
   // energy per unit mass
   double e_l = W_l.pressure/W_l.density/(gamma-1.0) + 0.5 * W_l.velocity.square();
   double e_r = W_r.pressure/W_r.density/(gamma-1.0) + 0.5 * W_r.velocity.square();
   
   // roe average
   double h = h_l * fact_l + h_r * fact_r;
   double c = sqrt( (gamma-1.0) * (h - 0.5*velocity.square()) );
   
   // speed of sound at l and r
   double s_l = min(vel_normal-c, v_l_normal-c_l);
   double s_r = max(vel_normal+c, v_r_normal+c_r);
   
   // speed of contact
   double s_m = (W_l.pressure - W_r.pressure
                 - W_l.density * v_l_normal * (s_l-v_l_normal)
                 + W_r.density * v_r_normal * (s_r-v_r_normal))
   /(W_r.density*(s_r-v_r_normal) - W_l.density*(s_l-v_l_normal));
   
   // Pressure at right and left (Pressure_j=Pressure_i) side of contact surface
   double pStar = W_r.density * (v_r_normal-s_r)*(v_r_normal-s_m) + W_r.pressure;
   
   if (s_m >= 0.0) {
      if (s_l > 0.0)
      {
         flux.mass_flux = W_l.density * v_l_normal;
         flux.momentum_flux = W_l.velocity * W_l.density * v_l_normal + unormal * W_l.pressure;
         flux.energy_flux = W_l.density * h_l * v_l_normal;
      }
      else
      {
         double invSLmSs = 1.0/(s_l-s_m);
         double sLmuL = s_l-v_l_normal;
         double rhoSL = W_l.density*sLmuL*invSLmSs;
         Vector rhouSL = (W_l.velocity*W_l.density*sLmuL+unormal*(pStar-W_l.pressure))*invSLmSs;
         double eSL = (sLmuL*e_l*W_l.density-W_l.pressure*v_l_normal+pStar*s_m)*invSLmSs;
         
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
         double sRmuR = s_r-v_r_normal;
         double rhoSR = W_r.density*sRmuR*invSRmSs;
         Vector rhouSR = (W_r.velocity*W_r.density*sRmuR+unormal*(pStar-W_r.pressure))*invSRmSs;
         double eSR = (sRmuR*e_r*W_r.density-W_r.pressure*v_r_normal+pStar*s_m)*invSRmSs;
         
         flux.mass_flux = rhoSR*s_m;
         flux.momentum_flux = rhouSR*s_m + unormal*pStar;
         flux.energy_flux = (eSR+pStar)*s_m;
      }
      else
      {
         flux.mass_flux = W_r.density*v_r_normal;
         flux.momentum_flux = W_r.velocity*W_r.density*v_r_normal + unormal * W_r.pressure;
         flux.energy_flux = W_r.density * h_r * v_r_normal;
      }
   }
   
   flux *= area;
}
