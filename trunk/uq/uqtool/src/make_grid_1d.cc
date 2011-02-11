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
//
// Linear element:
// Two linear elements and 3 dof
//   dofs    :          0-----1-----2
//   elements:          |--0--|--1--|
//
// Quadratic element:
// Two quadratic elements and 5 dof
//   dofs    :          0--1--2--3--4
//   elements:          |--0--|--1--|
template <>
void UQProblem<1>::make_grid ()
{
   cout << "Creating stochastic grid in 1-D ... ";
   
   // Number of initial samples
   unsigned int n_sample;
   if(order == 1)
      n_sample = 3;
   else if(order == 2)
      n_sample = 5;
   else
   {
      cout << "make_grid: Unknown initial order = " 
           << order << endl;
      abort ();
   }
   
   // Initial samples uniformly distributed
   double dx = (pdf_data.x_max[0] - pdf_data.x_min[0])/(n_sample-1);
   for(unsigned int i=0; i<n_sample; ++i)
   {
      Sample<1> new_sample (n_var, n_cell, n_moment, i);
      new_sample.x[0] = pdf_data.x_min[0] + i * dx;
      sample.push_back (new_sample);
   }
   
   // Initial elements, all second order
   unsigned int n_element = 2;
   for(unsigned int i=0; i<n_element; ++i)
   {      
      Element<1> new_element(order, n_moment, n_cell, i);
      grid.element.push_back (new_element);
   }
   
   // Map element dof to samples
   if(order==1)
   {
      // First element
      grid.element[0].idof[0] = 0; // left vertex
      grid.element[0].idof[1] = 1; // right vertex
      
      // Second element
      grid.element[1].idof[0] = 1; // left vertex
      grid.element[1].idof[1] = 2; // right vertex
   }
   else if(order==2)
   {
      // First element
      grid.element[0].idof[0] = 0; // left vertex
      grid.element[0].idof[1] = 2; // right vertex
      grid.element[0].idof[2] = 1; // middle
      
      // Second element
      grid.element[1].idof[0] = 2; // left vertex
      grid.element[1].idof[1] = 4; // right vertex
      grid.element[1].idof[2] = 3; // middle
   }
   
   grid.reinit_dof (sample);
   
   cout << "Done\n";
}
