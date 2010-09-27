#include "matrix.h"

using namespace std;

// empty constructor
Matrix::Matrix ()
{
   nrow = 0;
   ncol = 0;
}

// constructor
Matrix::Matrix (const unsigned int nrow, const unsigned int ncol)
   :
   nrow (nrow),
   ncol (ncol)
{
   data = new double[nrow*ncol];
}

// destructor
Matrix::~Matrix ()
{
   if(nrow*ncol > 0)
      delete [] data;
}

// access matrix element (i,j)
double& Matrix::operator() (const unsigned int i, const unsigned int j) const
{
   return data[i + nrow*j];
}

// access matrix element (i,j)
double& Matrix::operator() (const unsigned int i, const unsigned int j)
{
   return data[i + nrow*j];
}

void Matrix::allocate (const unsigned int ni, const unsigned int nj)
{
   nrow = ni;
   ncol = nj;
   data = new double[nrow*ncol];
}
