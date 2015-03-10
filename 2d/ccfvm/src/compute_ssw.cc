#include <vector>
#include "fv.h"

using namespace std;

// Compute shock indicator of Liou
void FiniteVolume::compute_ssw()
{
   // compute only if flux = kepes_roe since others do not use ssw
   // This is used to implement carbuncle fix.
   if(param.material.flux_scheme != Material::kepes_roe &&
      param.material.flux_scheme != Material::kepes_roe2) return;

   for(unsigned int i=0; i<grid.n_vertex; ++i)
      ssw[i] = 0;

   // Loop over interior faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      unsigned int cl = grid.face[i].vertex[0];
      unsigned int cr = grid.face[i].vertex[1];

      double ui = primitive[cl].velocity * grid.face[i].normal;
      double uj = primitive[cr].velocity * grid.face[i].normal;

      double ai = param.material.sound_speed(primitive[cl]);
      double aj = param.material.sound_speed(primitive[cr]);

      double area = grid.face[i].normal.norm();

      double l1i = ui - ai * area;
      double l1j = uj - aj * area;

      double l5i = ui + ai * area;
      double l5j = uj + aj * area;

      if(l1i>0 && l1j<0) ssw[cl] = ssw[cr] = 1;
      if(l5i>0 && l5j<0) ssw[cl] = ssw[cr] = 1;
   }
}
