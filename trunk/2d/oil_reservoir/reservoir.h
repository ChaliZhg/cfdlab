#ifndef __RESERVOIR_H__
#define __RESERVOIR_H__

#include "matrix.h"
#include "grid.h"

// numerical flux functions
double s_num_flux ();
double c_num_flux ();

// Class for reservoir problem
class ReservoirProblem
{
   public:
      ReservoirProblem () {};
      ~ReservoirProblem () {};
      void run ();

   private:
      double  final_time, dt;
      Grid    grid;
      Matrix  saturation;
      Matrix  concentration;
      Matrix  pressure;

      void make_grid ();
      void initialize ();
      void residual (Matrix&, Matrix&);
      void solve ();
      void output (const unsigned int);

};

#endif
