#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// Local Lax-Friedrichs flux function
//------------------------------------------------------------------------------
void Material::lxf_flux (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         Flux& flux) const
{
   
   double area = normal.norm();
   Vector unit_normal = normal / area;

   double left_density  = Density (left);
   double right_density = Density (right);

   // Enthalpy
   double E_left  = left.pressure/(gamma-1.0)
                  + 0.5 * left_density * left.velocity.square();
   double E_right = right.pressure/(gamma-1.0)
                  + 0.5 * right_density * right.velocity.square();

   double vel_left_normal  = left.velocity  * unit_normal;
   double vel_right_normal = right.velocity * unit_normal;

   // Left flux
   flux.mass_flux = left_density * vel_left_normal;
   flux.momentum_flux = unit_normal * left.pressure +
                        left.velocity * left_density * vel_left_normal;
   flux.energy_flux = (E_left + left.pressure) * vel_left_normal;

   // Right flux
   flux.mass_flux += right_density * vel_right_normal;
   flux.momentum_flux += unit_normal * right.pressure +
                        right.velocity * right_density * vel_right_normal;
   flux.energy_flux += (E_right + right.pressure) * vel_right_normal;

   double c_left = sound_speed (left);
   double c_right = sound_speed (right);

   double amax = max( fabs(vel_left_normal)  + c_left, 
                      fabs(vel_right_normal) + c_right );

   flux.mass_flux -= (right_density - left_density) * amax;
   flux.momentum_flux -= (right.velocity * right_density - left.velocity * left_density) * amax;
   flux.energy_flux -= (E_right - E_left) * amax;

   flux *= 0.5 * area;
}
