/*
 *  quadrature_1d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 16/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "quadrature.h"
#include "dunavant.h"

using namespace std;

// For given element, compute quadrature points and weights in 1-D
// Gauss quadrature rule
template <>
void Quadrature<2>::reinit (const Element<2>& element)
{  
   // Free old memory that may have been allocated
   free_memory ();
   
   // Quadrature must be atleast exact for polynomials of
   degree = 2 * element.order + 1;
   
   double* tri_xy = new double[2*3];
   tri_xy[0] = element.dof[0]->x[0];
   tri_xy[1] = element.dof[0]->x[1];
   tri_xy[2] = element.dof[1]->x[0];
   tri_xy[3] = element.dof[1]->x[1];
   tri_xy[4] = element.dof[2]->x[0];
   tri_xy[5] = element.dof[2]->x[1];
   
   n_point = dunavant_order_num ( degree );
   allocate_memory ();
   
   // Reference coordinates of quadrature nodes
   double* xytab = new double[2*n_point];
   
   // Physical Coordinates of quadrature nodes
   double* xytab_phy = new double[2*n_point];
   
   dunavant_rule ( degree, n_point, xytab, weight );
   reference_to_physical_t3 ( tri_xy, n_point, xytab, xytab_phy );

   double area = triangle_area (tri_xy);
   
   for(unsigned int i=0; i<n_point; ++i)
   {
      q_point[i][0] = xytab_phy[0+2*i];
      q_point[i][1] = xytab_phy[1+2*i];
      weight[i] *= area;
   }
}

