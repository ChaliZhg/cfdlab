/*
 *  write_sol.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 16/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>

using namespace std;

void write_sol (const char* filename, 
                const unsigned int n_var,
                const unsigned int n_cell,
                const double* data)
{
   ofstream fo;
   fo.open (filename);
   fo.precision (15);
   fo.setf (ios::scientific);
   
   unsigned int c = 0;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      for(unsigned int j=0; j<n_var; ++j)
         fo << data[c++] << " ";
      fo << endl;
   }
   fo.close ();
}