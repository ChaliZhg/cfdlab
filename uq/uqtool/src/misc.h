/*
 *  misc.h
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 07/06/11.
 *  Copyright 2011 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#ifndef __MISC_H__
#define __MISC_H__

#include <cmath>

// Compute L2 norm of x
inline
double norm (const std::valarray<double>& x)
{
   int n = x.size();
   
   double result = 0.0;
   
   for(int i=0; i<n; ++i)
      result += x[i] * x[i];
   
   return std::sqrt (result);
}

// Compute area of 2-d triangle
inline
double tri_area (const std::valarray<double>& x0,
                 const std::valarray<double>& x1,
                 const std::valarray<double>& x2)
{
   std::valarray<double> d01 = x1 - x0;
   std::valarray<double> d12 = x2 - x0;
   
   return 0.5 * std::fabs( d01[0] * d12[1] - d01[1] * d12[0] );
}

#endif