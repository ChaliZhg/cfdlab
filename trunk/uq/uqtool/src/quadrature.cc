/*
 *  quadrature.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */


#include <iostream>
#include <cassert>
#include "quadrature.h"

using namespace std;

// Constructor
template <int dim>
Quadrature<dim>::Quadrature ()
{
   n_point = 0;
}

// Destructor
template <int dim>
Quadrature<dim>::~Quadrature()
{
   free_memory ();
}

// Free memory for weights and quadrature points
template <int dim>
void Quadrature<dim>::free_memory ()
{
   if(n_point > 0)
   {
      // Deallocate memory for quadrature points
      for(unsigned int i=0; i<n_point; ++i)
         delete [] q_point[i];
      
      delete [] q_point;
      delete [] weight;
      
      n_point = 0;
   }
}

// Allocate memory for weights and quadrature points
template <int dim>
void Quadrature<dim>::allocate_memory ()
{
   assert (n_point > 0);
   
   weight  = new double  [n_point];
   q_point = new double* [n_point];
   for(unsigned int i=0; i<n_point; ++i)
      q_point[i] = new double [dim];
}

// To avoid linker errors
template class Quadrature<1>;
template class Quadrature<2>;