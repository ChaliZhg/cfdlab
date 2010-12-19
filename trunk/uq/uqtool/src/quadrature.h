/*
 *  quadrature.h
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#ifndef __QUADRATURE_H__
#define __QUADRATURE_H__

#include "grid.h"

template <int dim>
class Quadrature
{
   public:
      Quadrature ();
      ~Quadrature ();
      void reinit (const Element<dim>& element);
      void allocate_memory ();
      void free_memory ();
      
      unsigned int degree;  // polynomial degree of exactness
      unsigned int n_point; // no. of quadrature points
      double**     q_point; // location of quadrature points
      double*      weight;  // weight for each quadrature point
};

#endif