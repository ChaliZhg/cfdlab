/*
 *  pdf.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "pdf.h"
#include "prob.h"

using namespace std;

// Computes probability of dim independent random variables
template <int dim>
double PDFData<dim>::get_pdf (const double* x) const
{
   double pdf = 1.0, p, xx;
   
   for(unsigned int i=0; i<dim; ++i)
   {
      switch (type[i])
      {
         // Uniform random variable
         case (int)PDFType::uniform:
            if(x[i]>=x_min[i] && x[i]<=x_max[i])
               p = 1.0/(x_max[i] - x_min[i]);
            else
               p = 0.0;
            break;
            
         // Normal/gaussian random variable
         case (int)PDFType::normal:
            p = exp(-0.5 * (x[i] - mean[i]) * (x[i] - mean[i]) / variance[i]);
            p = p / sqrt(2.0 * M_PI * variance[i]);
            break;

         // Beta random variable
         case (int)PDFType::beta:
            xx = (x[i] - x_min[i])/(x_max[i] - x_min[i]);
            // Prevent xx going beyond [0,1] due to round-off
            xx = max(0.0, xx);
            xx = min(1.0, xx);
            p = Probability::beta_pdf (xx, alpha[i], beta[i]);
            break;

         case (int)PDFType::lognormal:
            cout << "lognormal not implemented\n";
            abort ();

         default:
            cout << "Unknown PDF !!!" << endl;
            abort ();
      }
      pdf *= p;
   }
   
   return pdf;
}

// To avoid linker error
template class PDFData<1>;
template class PDFData<2>;
