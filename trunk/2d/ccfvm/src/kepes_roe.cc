#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// KEPS flux with entropy dissipation
//------------------------------------------------------------------------------
void Material::kepes_roe_flux (const PrimVar& left,
                               const PrimVar& right,
                               const Vector& normal,
                               Flux& flux) const
{
   static const double BETA = 1.0/6.0;
   
   double area = normal.norm();
   Vector unit_normal = normal / area;

   double rhol = left.density;
   double rhor = right.density;
   double rho = logavg( rhol, rhor );
   double vel2= 0.5 * (left.velocity.square() + right.velocity.square());
   double betal = 0.5 * left.density / left.pressure;
   double betar = 0.5 * right.density / right.pressure;
   double beta = logavg(betal, betar);
   double a   = sqrt(0.5 * gamma / beta);

   double p     = 0.5 * (rhol + rhor) / (betal + betar);

   Vector vel = (left.velocity + right.velocity) / 2.0;
   double vel_normal = vel * unit_normal;

   // central flux
   flux.mass_flux = rho * vel_normal;
   flux.momentum_flux = unit_normal * p + vel * flux.mass_flux;
   flux.energy_flux = 0.5 * ( 1.0/((gamma-1.0)*beta) - vel2) * flux.mass_flux + 
                      flux.momentum_flux *  vel;

   // entropy dissipation
   // eigenvectors
   double H = a*a/(gamma-1.0) + 0.5*vel.square();
   double v1 = vel.x * unit_normal.y - vel.y * unit_normal.x;
   double R[][4] = {
      {            1,             1,               0,              1                      },
      {vel.x - a*unit_normal.x, vel.x,             unit_normal.y,  vel.x + a*unit_normal.x},
      {vel.y - a*unit_normal.y, vel.y,            -unit_normal.x,  vel.y + a*unit_normal.y},
      {H     - a*vel_normal,    0.5*vel.square(),  v1,             H     + a*vel_normal   }
   };

   // eigenvalues
   double vnl = left.velocity  * unit_normal;
   double vnr = right.velocity * unit_normal;
   double al  = sound_speed (left);
   double ar  = sound_speed (right);
   double LambdaL[] = { vnl - al, vnl, vnl, vnl + al };
   double LambdaR[] = { vnr - ar, vnr, vnr, vnr + ar };
   double l2, l3;
   l2 = l3 = fabs(vel_normal);
   double Lambda[]  = { fabs(vel_normal - a) + BETA*fabs(LambdaL[0]-LambdaR[0]), 
                        l2,
                        l3,
                        fabs(vel_normal + a) + BETA*fabs(LambdaL[3]-LambdaR[3])};

   double S[] = { 0.5*rho/gamma, (gamma-1.0)*rho/gamma, p, 0.5*rho/gamma };
   double D[] = { Lambda[0]*S[0], 
                  Lambda[1]*S[1],
                  Lambda[2]*S[2],
                  Lambda[3]*S[3]
                };

   // jump in entropy: s = log(p) - gamma*log(rho)
   double ds    = log(right.pressure/left.pressure) -
                  gamma * log(right.density/left.density);
   // Jump in entropy variables
   double dV[] = { -ds/(gamma-1.0) - 
                       (betar*right.velocity.square() - betal*left.velocity.square()),
                    2.0*(betar*right.velocity.x - betal*left.velocity.x),
                    2.0*(betar*right.velocity.y - betal*left.velocity.y),
                   -2.0*(betar - betal) };

   // DRT = D * R^T
   double DRT[4][4];
   for(unsigned int i=0; i<4; ++i)
      for(unsigned int j=0; j<4; ++j)
         DRT[i][j] = D[i]*R[j][i];

   // diffusive flux = R * Lambda * S * R^T * dV
   double Diff[] = {0.0, 0.0, 0.0, 0.0, 0.0};
   for(unsigned int i=0; i<4; ++i)
      for(unsigned int j=0; j<4; ++j)
         for(unsigned int k=0; k<4; ++k)
            Diff[i] += R[i][j] * DRT[j][k] * dV[k];

   flux.mass_flux       -= 0.5 * Diff[0];
   flux.momentum_flux.x -= 0.5 * Diff[1];
   flux.momentum_flux.y -= 0.5 * Diff[2];
   flux.energy_flux     -= 0.5 * Diff[3];

   flux *= area;
}
