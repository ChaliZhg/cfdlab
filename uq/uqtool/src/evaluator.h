/*
 *  evaluator.h
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 16/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__

#include "interpolate.h"

template <int dim>
class JREvaluator
{
public:
   JREvaluator (const unsigned int n_moment);
   ~JREvaluator ();
   void execute (const double* x,
                 const Interpolate<dim>& interpolate_formula);

   
   unsigned int n_moment;
   double* J;
   double* VdotR;
   double* RE;
};

#endif