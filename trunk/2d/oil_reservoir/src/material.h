#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <cmath>

#define harmonic_average(a,b)   (2.0*(a)*(b)/((a)+(b)))

const double viscosity_oil = 10.0;
const double density_water = 1.0;
const double density_oil   = 0.9;
extern double gravity;

// viscosity of water
double viscosity_water (const double& concentration);

// mobility of water
double mobility_water (const double& saturation, const double& concentration);

// mobility of oil
double mobility_oil (const double& saturation, const double& concentration);

// total mobility
double mobility_total (const double& saturation, const double& concentration);

// permeability of rock
double rock_permeability (const double& x, const double& y);


// viscosity of water as function of polymer concentration
inline
double viscosity_water (const double& concentration)
{
   return (1.0 + concentration);
}

// mobility of water
inline
double mobility_water (const double& saturation, const double& concentration)
{
   return saturation * saturation / viscosity_water (concentration);
}

// mobility of oil
inline
double mobility_oil (const double& saturation, const double& concentration)
{
   return (1.0 - saturation) * (1.0 - saturation) / viscosity_oil;
}

// total mobility
inline
double mobility_total (const double& saturation, const double& concentration)
{
   return mobility_water (saturation, concentration) +
          mobility_oil   (saturation, concentration);
}

// permeability of rock
inline
double rock_permeability (const double& x, const double& y)
{
   //return 1.0;
   return 1.0 + 0.5 * cos(6.0*M_PI*(x+0.2)) * cos(6.0*M_PI*y);
}

#endif
