/*
 *  pdf.h
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 15/12/10.
 *  Copyright 2010 TIFR-CAM, Bangalore. All rights reserved.
 *
 */

#ifndef __PDF_H__
#define __PDF_H__

#include <vector>
#include <string>

enum PDFType { uniform, normal };

template <int dim>
class PDFData
{
   public:
      PDFData () : x_name(dim) {};
      std::vector<std::string> x_name;
      double x_min[dim];
      double x_max[dim];
      PDFType type[dim];
      double  mean[dim];
   
      double get_pdf (const double* x) const;
};

#endif