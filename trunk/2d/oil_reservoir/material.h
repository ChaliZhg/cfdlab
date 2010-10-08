#define harmonic_average(a,b)   (2.0*(a)*(b)/((a)+(b)))

const double viscosity_oil = 10.0;

// viscosity of water
double viscosity_water (const double concentration);

// mobility of water
double mobility_water (const double saturation, const double concentration);

// mobility of oil
double mobility_oil (const double saturation, const double concentration);

// total mobility
double mobility_total (const double saturation, const double concentration);

// permeability of rock
double permeability (const double x, const double y);
