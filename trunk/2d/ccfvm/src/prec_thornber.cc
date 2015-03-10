#include <vector>
#include "fv.h"

using namespace std;

// Modify reconstructed states by Thornber scheme
void FiniteVolume::prec_thornber(vector<PrimVar>& state) const
{
   double ml = param.material.Mach (state[0]);
   double mr = param.material.Mach (state[1]);

   double z = min(1.0, max(ml, mr));

   Vector vl = state[0].velocity;
   Vector vr = state[1].velocity;

   state[0].velocity = (vl*(1.0+z) + vr*(1.0-z))/2.0;
   state[1].velocity = (vr*(1.0+z) + vl*(1.0-z))/2.0;
}
