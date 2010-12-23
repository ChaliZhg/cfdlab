/*
 *  interpolate.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include "interpolate.h"
#include "grid.h"

using namespace std;

// Constructor
template <int dim>
Interpolate<dim>::Interpolate (const unsigned int n_var,
                               const unsigned int n_cell)
   :
   n_var  (n_var),
   n_cell (n_cell)
{
   primal1  = new double [n_var * n_cell];
   adjoint1 = new double [n_var * n_cell];
   primal2  = new double [n_var * n_cell];
   adjoint2 = new double [n_var * n_cell];
   
   // Set pointer
   primal  = primal2;
   adjoint = adjoint2;
}

// Destructor
template <int dim>
Interpolate<dim>::~Interpolate ()
{
   primal  = NULL;
   adjoint = NULL;
   
   delete [] primal1;
   delete [] adjoint1;
   delete [] primal2;
   delete [] adjoint2;    
}

// Setup interpolation on given element
template <int dim>
void Interpolate<dim>::reinit (const Element<dim>& element)
{
   order = element.order;
   dof   = element.dof;
}

// To avoid linker errors
template class Interpolate<1>;