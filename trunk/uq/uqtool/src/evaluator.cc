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
JREvaluator<dim>::JREvaluator (const unsigned int n_moment)
   :
   n_moment (n_moment)
{
   J       = new double [n_moment];
   VdotR   = new double [n_moment];
   dVdotR  = new double [n_moment];
   dUdotAR = new double [n_moment];
}

// Destructor
template <int dim>
JREvaluator<dim>::~JREvaluator ()
{
   delete [] J;
   delete [] VdotR;
   delete [] dVdotR;
   delete [] dUdotAR;
}

template <int dim>
void JREvaluator<dim>::execute (const double* x,
                                const Interpolate<dim>& interpolate_formula)
{
   const unsigned int n_var  = interpolate_formula.n_var;
   const unsigned int n_cell = interpolate_formula.n_cell;
   const unsigned int n      = n_var * n_cell;
   
   char filename[64];
   ofstream fo;
   
   system ("mkdir -f EVAL");
   
   // Write random variables
   sprintf(filename, "EVAL/random.dat");
   fo.open (filename);
   for(unsigned int i=0; i<dim; ++i)
         fo << x[i] << endl;
   fo.close ();
   
   // Write primal solution
   sprintf(filename, "EVAL/primal.dat");
   write_sol (filename, n_var, n_cell, interpolate_formula.primal);

   // Write adjoint solution
   sprintf(filename, "EVAL/adjoint.dat");
   write_sol (filename, n_var, n_cell, interpolate_formula.adjoint);
   
   // Write primal2 - primal1
   double* dprimal = new double [n];
   for(unsigned int i=0; i<n; ++i)
      dprimal[i] = interpolate_formula.primal2[i] -
                   interpolate_formula.primal1[i];
   sprintf(filename, "EVAL/dprimal.dat");
   write_sol (filename, n_var, n_cell, dprimal);
   delete [] dprimal;
   
   // Write adjoint2 - adjoint1
   double* dadjoint = new double [n];
   for(unsigned int i=0; i<n; ++i)
      dadjoint[i] = interpolate_formula.adjoint2[i] -
                    interpolate_formula.adjoint1[i];
   sprintf(filename, "EVAL/dadjoint.dat");
   write_sol (filename, n_var, n_cell, dadjoint);
   delete [] dadjoint;
   
   // Run external simulator: DO NOT ITERATE
   system ("./runsolver.sh EVAL");
   
   ifstream fi;
   
   // Read objective functions
   fi.open ("EVAL/obj.dat");
   for(unsigned int i=0; i<n_moment; ++i)
      fi >> J[i] >> VdotR[i] >> dVdotR[i] >> dUdotAR[i];
   fi.close ();
   
}

// To avoid linker error
template class JREvaluator<1>;