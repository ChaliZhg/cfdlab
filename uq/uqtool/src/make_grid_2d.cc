/*
 *  make_grid_2d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 07/06/11.
 *  Copyright 2011 TIFR-CAM, Bangalore. All rights reserved.
 *
 */


#include <iostream>
#include "uq.h"

using namespace std;

// Make initial stochastic grid in 2-D
//
// Linear element:
// Two linear elements and 6 dof
//
//  1---------3
//  |        /|
//  |       / |
//  |      /  |
//  |     /   |
//  |    /    |
//  |   /     |
//  |  /      |
//  | /       |
//  |/        |
//  0---------2
//
// Quadratic element:
// Two quadratic elements and 9 dof
//
//  2----5----8
//  |        /|
//  |       / |
//  |      /  |
//  |     /   |
//  1    4    7
//  |   /     |
//  |  /      |
//  | /       |
//  |/        |
//  0----3----6
//
template <>
void UQProblem<2>::make_grid ()
{
   cout << "Creating stochastic grid in 2-D ... ";
   
   unsigned int nx;       // samples along each axis
   unsigned int n_sample; // total samples
   
   if(order == 1)
   {
      n_sample = 4;
      nx = 2;
   }
   else if(order == 2)
   {
      n_sample = 9;
      nx = 3;
   }
   else
   {
      cout << "make_grid: Unknown initial order = " 
           << order << endl;
      abort ();
   }
   
   // Initial samples uniformly distributed
   double dx[2];
   dx[0] = (pdf_data.x_max[0] - pdf_data.x_min[0])/(nx-1);
   dx[1] = (pdf_data.x_max[1] - pdf_data.x_min[1])/(nx-1);

   unsigned int n = 0;
   for(unsigned int i=0; i<nx; ++i)
      for(unsigned int j=0; j<nx; ++j)
      {
         Sample<2> new_sample (n_var, n_cell, n_moment, n);
         new_sample.x[0] = pdf_data.x_min[0] + i * dx[0];
         new_sample.x[1] = pdf_data.x_min[1] + j * dx[1];
         sample.push_back (new_sample);
         ++n;
      }
   
   // Initial elements, all second order
   unsigned int n_element = 2;
   for(unsigned int i=0; i<n_element; ++i)
   {      
      Element<2> new_element(order, n_moment, n_cell, i);
      grid.element.push_back (new_element);
   }
   
   // Map element dof to samples
   if(order==1)
   {
      // First element
      grid.element[0].idof[0] = 0;
      grid.element[0].idof[1] = 2;
      grid.element[0].idof[2] = 3;
      
      // Second element
      grid.element[1].idof[0] = 0;
      grid.element[1].idof[1] = 3;
      grid.element[1].idof[2] = 1;
   }
   else if(order==2)
   {
      // First element
      grid.element[0].idof[0] = 0;
      grid.element[0].idof[1] = 6;
      grid.element[0].idof[2] = 8;
      grid.element[0].idof[3] = 3;
      grid.element[0].idof[4] = 7;
      grid.element[0].idof[5] = 4;
      
      // Second element
      grid.element[1].idof[0] = 0;
      grid.element[1].idof[1] = 8;
      grid.element[1].idof[2] = 2;
      grid.element[1].idof[3] = 4;
      grid.element[1].idof[4] = 5;
      grid.element[1].idof[5] = 1;
   }
   
   grid.reinit_dof (sample);
   
   cout << "Done\n";
}
