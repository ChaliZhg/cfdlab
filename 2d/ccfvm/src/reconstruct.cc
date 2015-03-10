#include <cmath>
#include "fv.h"
#include "limiter.h"

#define EPSILON  1.0e-14
#define limit_albada(a,b)  max(0.0, (2*a*b + EPSILON)/(a*a + b*b + EPSILON))

using namespace std;

//------------------------------------------------------------------------------
// First order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_first
(
 const unsigned int& f,
 vector<PrimVar>&    state
) const
{
   // Left state
   unsigned int cl = grid.face[f].vertex[0];
   state[0] = primitive[cl];

   unsigned int cr = grid.face[f].vertex[1];
   state[1] = primitive[cr];
}

//------------------------------------------------------------------------------
// Second order Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_second
(
 const unsigned int& f,
 vector<PrimVar>&    state
) const
{
   unsigned int cl = grid.face[f].vertex[0];
   unsigned int cr = grid.face[f].vertex[1];

   Vector  dr    = grid.vertex[cr].coord - grid.vertex[cl].coord;
   PrimVar dprim = primitive[cr] - primitive[cl];

   // left state
   PrimVar dpriml;
   dpriml.temperature= 2.0 * (dT[cl] * dr) - dprim.temperature;
   dpriml.velocity.x = 2.0 * (dU[cl] * dr) - dprim.velocity.x;
   dpriml.velocity.y = 2.0 * (dV[cl] * dr) - dprim.velocity.y;
   dpriml.velocity.z = 2.0 * (dW[cl] * dr) - dprim.velocity.z;
   dpriml.pressure   = 2.0 * (dP[cl] * dr) - dprim.pressure;

   state[0] = primitive[cl] + (dpriml * (1.0-KKK) + dprim * (1.0+KKK)) * 0.25;

   // right state
   PrimVar dprimr;
   dprimr.temperature= 2.0 * (dT[cr] * dr) - dprim.temperature;
   dprimr.velocity.x = 2.0 * (dU[cr] * dr) - dprim.velocity.x;
   dprimr.velocity.y = 2.0 * (dV[cr] * dr) - dprim.velocity.y;
   dprimr.velocity.z = 2.0 * (dW[cr] * dr) - dprim.velocity.z;
   dprimr.pressure   = 2.0 * (dP[cr] * dr) - dprim.pressure;

   state[1] = primitive[cr] - (dprimr * (1.0-KKK) + dprim * (1.0+KKK)) * 0.25;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_limited
(
 const unsigned int& f,
 vector<PrimVar>&    state
) const
{
   unsigned int cl = grid.face[f].vertex[0];
   unsigned int cr = grid.face[f].vertex[1];

   Vector  dr    = grid.vertex[cr].coord - grid.vertex[cl].coord;
   PrimVar dprim = primitive[cr] - primitive[cl];

   // left state
   PrimVar dpriml;
   dpriml.temperature= 2.0 * (dT[cl] * dr) - dprim.temperature;
   dpriml.velocity.x = 2.0 * (dU[cl] * dr) - dprim.velocity.x;
   dpriml.velocity.y = 2.0 * (dV[cl] * dr) - dprim.velocity.y;
   dpriml.velocity.z = 2.0 * (dW[cl] * dr) - dprim.velocity.z;
   dpriml.pressure   = 2.0 * (dP[cl] * dr) - dprim.pressure;

   PrimVar si = limited_slope(dpriml, dprim);
   state[0] = primitive[cl] + si;

   // right state
   PrimVar dprimr;
   dprimr.temperature= 2.0 * (dT[cr] * dr) - dprim.temperature;
   dprimr.velocity.x = 2.0 * (dU[cr] * dr) - dprim.velocity.x;
   dprimr.velocity.y = 2.0 * (dV[cr] * dr) - dprim.velocity.y;
   dprimr.velocity.z = 2.0 * (dW[cr] * dr) - dprim.velocity.z;
   dprimr.pressure   = 2.0 * (dP[cr] * dr) - dprim.pressure;

   PrimVar sj = limited_slope(dprimr, dprim);
   state[1] = primitive[cr] - sj;
}

//------------------------------------------------------------------------------
// Computed limited slope
//------------------------------------------------------------------------------
PrimVar FiniteVolume::limited_slope (const PrimVar& ul, const PrimVar& ur) const
{
   PrimVar s, result;

   s.temperature = limit_albada (ul.temperature, ur.temperature);
   s.velocity.x  = limit_albada (ul.velocity.x , ur.velocity.x);
   s.velocity.y  = limit_albada (ul.velocity.y , ur.velocity.y);
   s.velocity.z  = limit_albada (ul.velocity.z , ur.velocity.z);
   s.pressure    = limit_albada (ul.pressure   , ur.pressure);

   result.temperature = 0.25 * s.temperature *
                        ( (1.0-KKK*s.temperature)*ul.temperature + 
                          (1.0+KKK*s.temperature)*ur.temperature );

   result.velocity.x = 0.25 * s.velocity.x *
                        ( (1.0-KKK*s.velocity.x)*ul.velocity.x + 
                          (1.0+KKK*s.velocity.x)*ur.velocity.x );

   result.velocity.y = 0.25 * s.velocity.y *
                        ( (1.0-KKK*s.velocity.y)*ul.velocity.y + 
                          (1.0+KKK*s.velocity.y)*ur.velocity.y );

   result.velocity.z = 0.25 * s.velocity.z *
                        ( (1.0-KKK*s.velocity.z)*ul.velocity.z + 
                          (1.0+KKK*s.velocity.z)*ur.velocity.z );

   result.pressure = 0.25 * s.pressure *
                        ( (1.0-KKK*s.pressure)*ul.pressure + 
                          (1.0+KKK*s.pressure)*ur.pressure );

   return result;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_minmod
(
 const unsigned int& f,
 vector<PrimVar>&    state
 ) const
{
   unsigned int cl = grid.face[f].vertex[0];
   unsigned int cr = grid.face[f].vertex[1];
   
   Vector  dr    = grid.vertex[cr].coord - grid.vertex[cl].coord;
   PrimVar dprim = primitive[cr] - primitive[cl];
   
   // left state
   PrimVar dpriml;
   dpriml.temperature= 2.0 * (dT[cl] * dr) - dprim.temperature;
   dpriml.velocity.x = 2.0 * (dU[cl] * dr) - dprim.velocity.x;
   dpriml.velocity.y = 2.0 * (dV[cl] * dr) - dprim.velocity.y;
   dpriml.velocity.z = 2.0 * (dW[cl] * dr) - dprim.velocity.z;
   dpriml.pressure   = 2.0 * (dP[cl] * dr) - dprim.pressure;
   
   PrimVar si = minmod_slope(dpriml, dprim);
   state[0] = primitive[cl] + si;
   
   // right state
   PrimVar dprimr;
   dprimr.temperature= 2.0 * (dT[cr] * dr) - dprim.temperature;
   dprimr.velocity.x = 2.0 * (dU[cr] * dr) - dprim.velocity.x;
   dprimr.velocity.y = 2.0 * (dV[cr] * dr) - dprim.velocity.y;
   dprimr.velocity.z = 2.0 * (dW[cr] * dr) - dprim.velocity.z;
   dprimr.pressure   = 2.0 * (dP[cr] * dr) - dprim.pressure;
   
   PrimVar sj = minmod_slope(dprimr, dprim);
   state[1] = primitive[cr] - sj;
}

//------------------------------------------------------------------------------
// Computed limited slope
//------------------------------------------------------------------------------
PrimVar FiniteVolume::minmod_slope (const PrimVar& ul, const PrimVar& ur) const
{
   PrimVar s;
      
   s.temperature = 0.5 * minmod (ul.temperature, ur.temperature);
   s.velocity.x  = 0.5 * minmod (ul.velocity.x , ur.velocity.x);
   s.velocity.y  = 0.5 * minmod (ul.velocity.y , ur.velocity.y);
   s.velocity.z  = 0.5 * minmod (ul.velocity.z , ur.velocity.z);
   s.pressure    = 0.5 * minmod (ul.pressure   , ur.pressure);
   
   return s;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct_minmax
(
 const unsigned int& f,
 vector<PrimVar>&    state
 ) const
{
   unsigned int cl = grid.face[f].vertex[0];
   unsigned int cr = grid.face[f].vertex[1];
   
   Vector  dr    = grid.vertex[cr].coord - grid.vertex[cl].coord;
   
   // left state
   state[0].temperature = primitive[cl].temperature+ 0.5 * phi[cl].temperature * (dT[cl] * dr);
   state[0].velocity.x  = primitive[cl].velocity.x + 0.5 * phi[cl].velocity.x  * (dU[cl] * dr);
   state[0].velocity.y  = primitive[cl].velocity.y + 0.5 * phi[cl].velocity.y  * (dV[cl] * dr);
   state[0].velocity.z  = primitive[cl].velocity.z + 0.5 * phi[cl].velocity.z  * (dW[cl] * dr);
   state[0].pressure    = primitive[cl].pressure   + 0.5 * phi[cl].pressure    * (dP[cl] * dr);
   
   // right state
   state[1].temperature = primitive[cr].temperature- 0.5 * phi[cr].temperature * (dT[cr] * dr);
   state[1].velocity.x  = primitive[cr].velocity.x - 0.5 * phi[cr].velocity.x  * (dU[cr] * dr);
   state[1].velocity.y  = primitive[cr].velocity.y - 0.5 * phi[cr].velocity.x  * (dV[cr] * dr);
   state[1].velocity.z  = primitive[cr].velocity.z - 0.5 * phi[cr].velocity.x  * (dW[cr] * dr);
   state[1].pressure    = primitive[cr].pressure   - 0.5 * phi[cr].pressure    * (dP[cr] * dr);
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void FiniteVolume::reconstruct (const unsigned int& f,
                                vector<PrimVar>&    state) const
{
   switch(param.reconstruct_scheme)
   {
      // First order
      case Parameter::first:
         reconstruct_first (f, state);
         break;

      // Second order, MUSCL
      case Parameter::second:
         reconstruct_second (f, state);
         break;

      // MUSCL with van Albada limiter
      case Parameter::limited:
         reconstruct_limited (f, state);
         break;
         
      // MUSCL with minmod limiter
      case Parameter::minmod:
         reconstruct_minmod (f, state);
         break;

      // Barth-Jespersen or MinMax scheme
      case Parameter::bj:
      case Parameter::minmax:
         reconstruct_minmax (f, state);
         break;

      default:
         cout << "reconstruct: unknown reconstruction scheme = " 
              << param.reconstruct_scheme << endl;
         exit (0);
   }

   // Modify velocity by Thornber scheme for low mach case
   //prec_thornber(state);
}
