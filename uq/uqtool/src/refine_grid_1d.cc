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
template <>
void UQProblem<1>::refine_grid ()
{
   for(unsigned int i=0; i<grid.element.size(); ++i)
      if(grid.element[i].active && grid.element[i].refine_flag)
      {
         // TBD refine this element into two
      }
}
