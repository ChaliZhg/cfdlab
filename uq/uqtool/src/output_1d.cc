/*
 *  output_1d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 19/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <fstream>
#include "uq.h"

using namespace std;

// Save stochastic grid and functionals into file for
// visualization. Note: Spatial order may not be maintained
// need to sort based on x[0]
template <>
void UQProblem<1>::output (const unsigned int iter) const
{
   char filename[64];
   
   sprintf(filename, "RESULT/grid_%d.dat", iter);
   ofstream fo;
   fo.open (filename);
   fo.precision (15);
   
   for(unsigned int i=0; i<sample.size(); ++i)
   {
      fo << sample[i].x[0] << " ";
      for(unsigned int j=0; j<n_moment; ++j)
         fo << sample[i].J[j] << " ";
      fo << endl;
   }
   
   fo.close ();
   
   sprintf(filename, "RESULT/cc_%d.dat", iter);
   fo.open (filename);
   fo.precision (15);
   
   for(unsigned int i=0; i<grid.element.size(); ++i)
      if(grid.element[i].active)
      {
         double x = 0.5 * (grid.element[i].dof[0]->x[0] +
                           grid.element[i].dof[1]->x[0]);
         fo << x << " ";
         for(unsigned int j=0; j<n_moment; ++j)
            fo << grid.element[i].adj_cor[j] << " ";
         fo << endl;
      }
   
   fo.close ();
}
