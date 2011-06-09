/*
 *  interpolate.h
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

#include "grid.h"

template <int dim>
class Interpolate
{
public:
   Interpolate (const unsigned int n_var,
                const unsigned int n_cell);
   ~Interpolate ();
   void reinit (const Element<dim>& element);
   void execute (const double* x);
   double lagrange (const unsigned int iorder,
                    const unsigned int vertex,
                    const double*      x);
   
   unsigned int order;
   unsigned int n_var;
   unsigned int n_cell;
   double* primal1;  // low order interpolation
   double* adjoint1; // low order interpolation
   double* primal2;  // high order interpolation
   double* adjoint2; // high order interpolation
   double* primal;
   double* adjoint;
   std::vector<typename Sample<dim>::Sample*> dof; 

};

#endif