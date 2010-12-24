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
}
