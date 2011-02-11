/*
 *  refine_grid_1d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 17/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include "uq.h"

using namespace std;

// Refine the stochastic grid in 1-D
// * = new samples
//
// Linear element:
// old element:  0-----------1
// new element:  0-----*-----1
//
// Quadratic element:
// old element:  0-----2-----1
// new element:  0--*--2--*--1
template <>
void UQProblem<1>::refine_grid ()
{
   unsigned int n_element = grid.element.size();
   
   for(unsigned int i=0; i<n_element; ++i)
      if(grid.element[i].refine_flag)
      {         
         // Refine this element into two
         cout << "Refining element = " << i << endl;
          
         if(grid.element[i].order == 1)
         {
            unsigned int s0 = grid.element[i].idof[0];
            unsigned int s1 = grid.element[i].idof[1];
            
            // One new sample
            Sample<1> new_sample1 (n_var, n_cell, n_moment, sample.size());
            new_sample1.x[0] = (sample[s0].x[0] + sample[s1].x[0]) / 2.0;
            sample.push_back (new_sample1);
            
            // Two new linear elements
            Element<1> new_element1 (1, n_moment, n_cell, grid.element.size());
            new_element1.idof[0] = s0;
            new_element1.idof[1] = sample.size() - 1;
            grid.element.push_back (new_element1);
            
            Element<1> new_element2 (1, n_moment, n_cell, grid.element.size());
            new_element2.idof[0] = sample.size() - 1;
            new_element2.idof[1] = s1;
            grid.element.push_back (new_element2);
         }
         else if(grid.element[i].order == 2)
         {
            unsigned int s0 = grid.element[i].idof[0];
            unsigned int s1 = grid.element[i].idof[1];
            unsigned int s2 = grid.element[i].idof[2];
            
            // Two new samples
            Sample<1> new_sample1 (n_var, n_cell, n_moment, sample.size());
            new_sample1.x[0] = (sample[s0].x[0] + sample[s2].x[0]) / 2.0;
            sample.push_back (new_sample1);

            Sample<1> new_sample2 (n_var, n_cell, n_moment, sample.size());
            new_sample2.x[0] = (sample[s1].x[0] + sample[s2].x[0]) / 2.0;
            sample.push_back (new_sample2);
            
            // Two new quadratic elements
            Element<1> new_element1 (2, n_moment, n_cell, grid.element.size());
            new_element1.idof[0] = s0;
            new_element1.idof[1] = s2;
            new_element1.idof[2] = sample.size() - 2;
            grid.element.push_back (new_element1);
            
            Element<1> new_element2 (2, n_moment, n_cell, grid.element.size());
            new_element2.idof[0] = s2;
            new_element2.idof[1] = s1;
            new_element2.idof[2] = sample.size() - 1;
            grid.element.push_back (new_element2);
         }
         else
         {
            cout << "Cannot refine grid for order =" 
                 << grid.element[i].order
                 << endl;
            abort ();
         }

         grid.element[i].active = false;
         grid.element[i].refine_flag = false;
      }   
}
