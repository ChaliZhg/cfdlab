#ifndef __RESERVOIR_H__
#define __RESERVOIR_H__

#include "matrix.h"
#include "grid.h"

// numerical flux functions
double c_num_flux ();

// Class for reservoir problem
class ReservoirProblem
{
   public:
      ReservoirProblem () {};
      ~ReservoirProblem () {};
      void run ();

   private:
      unsigned int max_iter;
      double  cfl, final_time, dt;
      double  max_velocity;
      Grid    grid;
      Matrix  saturation;
      Matrix  concentration;
      Matrix  pressure;

      void make_grid ();
      void initialize ();
      void residual (Matrix&, Matrix&);
      void solve ();
      void output (const unsigned int);

      std::vector<double> num_flux 
         (const unsigned int, const unsigned int,
          const unsigned int, const unsigned int);
      void updateGhostCells ();
      void findMinMax ();

};

#endif
