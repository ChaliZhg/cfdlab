/*
 *  refine_grid_1d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 17/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include "uq.h"

// Refine the stochastic grid in 1-D
// TBD We are assuming second order elements
// * = new samples
// old element:  0-----2-----1
// new element:  0--*--2--*--1
template <>
void UQProblem<1>::refine_grid ()
{
   for(unsigned int i=0; i<grid.element.size(); ++i)
      if(grid.element[i].active && grid.element[i].refine_flag)
      {
         // Refine this element into two
         
         // Two new samples
         Sample<1> new_sample1 (n_var, n_cell, n_moment, sample.size());
         new_sample1.x[0] = (grid.element[i].dof[0]->x[0] + 
                             grid.element[i].dof[2]->x[0]) / 2.0;
         sample.push_back (new_sample1);

         Sample<1> new_sample2 (n_var, n_cell, n_moment, sample.size());
         new_sample2.x[0] = (grid.element[i].dof[1]->x[0] + 
                             grid.element[i].dof[2]->x[0]) / 2.0;
         sample.push_back (new_sample2);
         
         // Two new elements
         Element<1> new_element1 (2, n_moment);
         new_element1.parent = &grid.element[i];
         new_element1.dof[0] = grid.element[i].dof[0];
         new_element1.dof[1] = grid.element[i].dof[2];
         new_element1.dof[2] = &sample[sample.size()-2];
         grid.element.push_back (new_element1);

         Element<1> new_element2 (2, n_moment);
         new_element2.parent = &grid.element[i];
         new_element2.dof[0] = grid.element[i].dof[2];
         new_element2.dof[1] = grid.element[i].dof[1];
         new_element2.dof[2] = &sample[sample.size()-1];
         grid.element.push_back (new_element2);
         
         grid.element[i].active = false;
         grid.element[i].refine_flag = false;
      }
}
