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

#include <valarray>
#include "interpolate.h"

template <int dim>
class JREvaluator
{
public:
   JREvaluator (const std::vector<std::string> x_name,
                const unsigned int n_moment,
                const unsigned int n_cell,
                char* template_dir_);
   ~JREvaluator ();
   void execute (const double* x,
                 const Interpolate<dim>& interpolate_formula);

   std::vector<std::string> x_name;
   unsigned int n_moment;
   unsigned int n_cell;
   char template_dir[64];
   std::valarray<double> J;
   std::valarray<double> VdotR;
   std::valarray<double> RE;
   std::valarray<double> RE_array;
};

#endif