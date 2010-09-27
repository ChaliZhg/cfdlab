#ifndef __MATRIX_H__
#define __MATRIX_H__

class Matrix
{
   public:
      Matrix ();
      Matrix (const unsigned int nrow, const unsigned ncol);
      ~Matrix ();
      double& operator() (const unsigned int i, const unsigned int j) const;
      double& operator() (const unsigned int i, const unsigned int j);
      void allocate (const unsigned int, const unsigned int);

   private:
      unsigned int nrow, ncol;
      double* data;
};

#endif
