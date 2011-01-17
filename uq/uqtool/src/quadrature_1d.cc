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

using namespace std;

// For given element, compute quadrature points and weights in 1-D
// Gauss quadrature rule
template <>
void Quadrature<1>::reinit (const Element<1>& element)
{  
   // Free old memory that may have been allocated
   free_memory ();
   
   // Quadrature must be atleast exact for polynomials of
   // degree = 2 * element.order + 1;
      
   // Set number of quadrature points
   if(element.order == 1)
   {
      degree  = 3;
      n_point = 2;
      allocate_memory ();
      q_point[0][0] = -1.0/sqrt(3.0); weight[0] = 1.0;
      q_point[1][0] = +1.0/sqrt(3.0); weight[1] = 1.0;
   }
   else if(element.order == 2)
   {
      degree  = 7;
      n_point = 4;
      allocate_memory ();
      double q1 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
      double q2 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
      double w1 = (18.0 + sqrt(30.0))/36.0;
      double w2 = (18.0 - sqrt(30.0))/36.0;

      q_point[0][0] = -q2; weight[0] = w2;
      q_point[1][0] = -q1; weight[1] = w1;
      q_point[2][0] = +q1; weight[2] = w1;
      q_point[3][0] = +q2; weight[3] = w2;
   }
   else 
   {
      cout << "Quadrature<1> not defined for degree =" << degree
           << endl;
      abort ();
   }
      
   // Length of element
   double length = element.dof[1]->x[0] - element.dof[0]->x[0];
   
   // Rule is given for [-1,+1]. Hence divide weight by 2
   // Convert q_point from [-1,+1] to actual element interval
   for(unsigned int i=0; i<n_point; ++i)
   {
      q_point[i][0] = element.dof[0]->x[0] + 
                      0.5 * length * ( q_point[i][0] + 1.0 );
      weight[i] *= length / 2.0;
   }
}

