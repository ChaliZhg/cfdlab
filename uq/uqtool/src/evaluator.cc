/*
 *  evaluator.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 16/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include "evaluator.h"

using namespace std;

// Function declaration
void write_sol (const char* filename, 
                const unsigned int n_var,
                const unsigned int n_cell,
                const double* data);

// Constructor
template <int dim>
JREvaluator<dim>::JREvaluator (const vector<string> x_name,
                               const unsigned int n_moment,
                               const unsigned int n_cell)
   :
   x_name   (x_name),
   n_moment (n_moment),
   n_cell   (n_cell)
{
   J.resize        (n_moment);
   VdotR.resize    (n_moment);
   RE.resize       (n_moment);
   RE_array.resize (n_moment * n_cell);
}

// Destructor
template <int dim>
JREvaluator<dim>::~JREvaluator ()
{
}

template <int dim>
void JREvaluator<dim>::execute (const double* x,
                                const Interpolate<dim>& interpolate_formula)
{
   const unsigned int n_var  = interpolate_formula.n_var;
   const unsigned int n      = n_var * n_cell;
   assert (n_cell == interpolate_formula.n_cell);

   char filename[64];
   ofstream fo;
   
   system ("mkdir -p RESULT/EVAL");
   
   // Write random variables
   sprintf(filename, "RESULT/EVAL/random.dat");
   fo.open (filename);
   fo.precision (15);
   fo.setf (ios::scientific);
   for(unsigned int i=0; i<dim; ++i)
         fo << "{" << x_name[i] << "}  " << x[i] << endl;
   fo.close ();
   
   // Write primal solution
   sprintf(filename, "RESULT/EVAL/primal.dat");
   write_sol (filename, n_var, n_cell, interpolate_formula.primal);

   // Write adjoint solution
   sprintf(filename, "RESULT/EVAL/adjoint.dat");
   write_sol (filename, n_var, n_cell, interpolate_formula.adjoint);
   
   // Write primal2 - primal1
   double* dprimal = new double [n];
   for(unsigned int i=0; i<n; ++i)
      dprimal[i] = interpolate_formula.primal2[i] -
                   interpolate_formula.primal1[i];
   sprintf(filename, "RESULT/EVAL/dprimal.dat");
   write_sol (filename, n_var, n_cell, dprimal);
   delete [] dprimal;
   
   // TBD adjoint will be different for each moment
   // Write adjoint2 - adjoint1
   double* dadjoint = new double [n];
   for(unsigned int i=0; i<n; ++i)
      dadjoint[i] = interpolate_formula.adjoint2[i] -
                    interpolate_formula.adjoint1[i];
   sprintf(filename, "RESULT/EVAL/dadjoint.dat");
   write_sol (filename, n_var, n_cell, dadjoint);
   delete [] dadjoint;
   
   // Run external simulator: DO NOT ITERATE
   system ("./runsolver.sh 2 RESULT/EVAL");
   
   ifstream fi;
   
   // Read objective functions
   fi.open ("RESULT/EVAL/obj.dat");
   for(unsigned int i=0; i<n_moment; ++i)
      fi >> J[i];
   fi.close ();
   
   double VdotR_array;

   unsigned int c = 0;
   for(unsigned int i=0; i<n_moment; ++i)
   {      
      VdotR[i] = 0.0;
      RE[i]    = 0.0;
      
      // TBD Filename will be different for each moment
      fi.open ("RESULT/EVAL/VdotR.dat");
      assert (fi.is_open());
      for(unsigned int j=0; j<n_cell; ++j)
      {
         fi >> VdotR_array;
         VdotR[i] += VdotR_array;
      }
      fi.close ();
      
      // TBD Filename will be different for each moment
      fi.open ("RESULT/EVAL/RE.dat");
      assert (fi.is_open());
      for(unsigned int j=0; j<n_cell; ++j)
      {
         fi >> RE_array[c];
         RE[i] += RE_array[c];
         ++c;
      }
      fi.close ();
   }
   
}

// To avoid linker error
template class JREvaluator<1>;