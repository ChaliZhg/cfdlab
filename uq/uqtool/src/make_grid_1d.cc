/*
 *  make_grid_1d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 16/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include "uq.h"

using namespace std;

// Make initial stochastic grid in 1-D
// Initial grid has two quadratic elements and 5 dof
//   dofs    :          0--1--2--3--4
//   elements:          |--0--|--1--|
template <>
void UQProblem<1>::make_grid ()
{
   cout << "Creating stochastic grid in 1-D ... ";
   
   // Initial 5 samples uniformly distributed
   n_sample = 5;
   double dx = (pdf_data.x_max[0] - pdf_data.x_min[0])/(n_sample-1);
   for(unsigned int i=0; i<n_sample; ++i)
   {
      Sample<1> new_sample (n_var, n_cell, n_moment, i);
      new_sample.x[0] = pdf_data.x_min[0] + i * dx;
      sample.push_back (new_sample);
   }
   
   // Initial elements, all second order
   grid.n_element = 2;
   for(unsigned int i=0; i<grid.n_element; ++i)
   {      
      Element<1> new_element(2, n_moment);
      grid.element.push_back (new_element);
   }
   
   // Map element dof to samples
   // First element
   grid.element[0].dof[0] = &sample[0]; // left vertex
   grid.element[0].dof[1] = &sample[2]; // right vertex
   grid.element[0].dof[2] = &sample[1]; // middle
   
   // Second element
   grid.element[1].dof[0] = &sample[2]; // left vertex
   grid.element[1].dof[1] = &sample[4]; // right vertex
   grid.element[1].dof[2] = &sample[3]; // middle
   
   cout << "Done\n";
}
