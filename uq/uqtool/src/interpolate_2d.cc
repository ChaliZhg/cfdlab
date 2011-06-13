/*
 *  interpolate_2d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 07/06/11.
 *  Copyright 2011 TIFR-CAM, Bangalore. All rights reserved.
 *
 */


#include <iostream>
#include <cassert>
#include "interpolate.h"
#include "misc.h"

using namespace std;

// Compute Lagrange basis function for triangle
template<>
double Interpolate<2>::lagrange (const unsigned int iorder,
                                 const unsigned int vertex,
                                 const double*      x)
{
   assert (iorder==1 || iorder==2);
   assert (vertex < 6);
   
   double area = tri_area (dof[0]->x, dof[1]->x, dof[2]->x);
   
   valarray<double> xp(2);
   xp[0] = x[0];
   xp[1] = x[1];
   
   double L[3];
   L[0] = tri_area (xp, dof[1]->x, dof[2]->x) / area;
   L[1] = tri_area (xp, dof[2]->x, dof[0]->x) / area;
   L[2] = tri_area (xp, dof[0]->x, dof[1]->x) / area;

   double result = 0.0;
   
   if(iorder == 1)
   {
      result = L[vertex];
   }
   else if (iorder == 2)
   {
      if(vertex<3)
         result = (2.0*L[vertex] - 1.0) * L[vertex];
      else if(vertex == 3)
         result = 4.0*L[0]*L[1];
      else if(vertex == 4)
         result = 4.0*L[1]*L[2];
      else if(vertex == 5)
         result = 4.0*L[2]*L[0];
   }
   
   return result;
}

// Lagrange interpolate of primal and adjoint solution in 1-D
template <>
void Interpolate<2>::execute (const double* x)
{
   assert (order == 1 || order == 2);
   
   double a[6];
   const unsigned int n = n_var * n_cell;
   
   // Linear interpolation
   for(unsigned int i=0; i<3; ++i)
      a[i] = lagrange(1, i, x);
   
   if(order == 1)
   {
      for(unsigned int i=0; i<n; ++i)
      {
         // low order is averaging
         primal1[i]  = ( dof[0]->primal[i] +
                         dof[1]->primal[i] +
                         dof[2]->primal[i] ) / 3.0;
         adjoint1[i] = ( dof[0]->adjoint[i] +
                         dof[1]->adjoint[i] +
                         dof[2]->adjoint[i] ) / 3.0;
         
         // high order is linear interpolation
         primal2[i]  = a[0] * dof[0]->primal[i] +
                       a[1] * dof[1]->primal[i] +
                       a[2] * dof[2]->primal[i];
         adjoint2[i] = a[0] * dof[0]->adjoint[i] +
                       a[1] * dof[1]->adjoint[i] +
                       a[2] * dof[2]->adjoint[i];
      }
   }
   else
   {
      // low order is linear interpolation
      for(unsigned int i=0; i<n; ++i)
      {
         primal1[i]  = a[0] * dof[0]->primal[i] +
                       a[1] * dof[1]->primal[i] +
                       a[2] * dof[2]->primal[i];
         adjoint1[i] = a[0] * dof[0]->adjoint[i] +
                       a[1] * dof[1]->adjoint[i] +
                       a[2] * dof[2]->adjoint[i];
      }
      
      // high order is Quadratic interpolation
      for(unsigned int i=0; i<6; ++i)
         a[i] = lagrange (2, i, x);
      
      for(unsigned int i=0; i<n; ++i)
      {
         primal2[i] = adjoint2[i] = 0.0;
         
         for(unsigned int j=0; j<6; ++j)
         {
            primal2[i]  += a[j] * dof[j]->primal[i];
            adjoint2[i] += a[j] * dof[j]->adjoint[i];
         }
      }
   }
   
}
