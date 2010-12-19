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
   
   // Quadrature must be exact for polynomials of "degree"
   degree = 2 * element.order + 1;
   
   // Set number of quadrature points
   if(degree == 3)
   {
      n_point = 2;
      allocate_memory ();
      q_point[0][0] = -1.0/sqrt(3.0); weight[0] = 1.0;
      q_point[1][0] = +1.0/sqrt(3.0); weight[1] = 1.0;
   }
   else if(degree == 5)
   {
      n_point = 3;
      allocate_memory ();
      q_point[0][0] = -sqrt(15.0)/5.0; weight[0] = 5.0/9.0;
      q_point[1][0] =  0.0;            weight[1] = 8.0/9.0;
      q_point[2][0] = +sqrt(15.0)/5.0; weight[2] = 5.0/9.0;
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
      q_point[i][0] = 
         2.0 * (q_point[i][0] - element.dof[0]->x[0]) / length - 1.0;
      weight[i] *= length / 2.0;
   }
}

