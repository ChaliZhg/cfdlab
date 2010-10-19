#include "material.h"

// viscosity of water as function of polymer concentration
double viscosity_water (const double concentration)
{
   return 1.0;
}

// mobility of water
double mobility_water (const double saturation, const double concentration)
{
   return saturation * saturation / viscosity_water (concentration);
}

// mobility of oil
double mobility_oil (const double saturation, const double concentration)
{
   return (1.0 - saturation) * (1.0 - saturation) / viscosity_oil;
}

// total mobility
double mobility_total (const double saturation, const double concentration)
{
   return mobility_water (saturation, concentration) +
          mobility_oil   (saturation, concentration);
}

// permeability of rock
double permeability (const double x, const double y)
{
   return 1.0;
}
