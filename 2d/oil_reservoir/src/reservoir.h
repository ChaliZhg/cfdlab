#ifndef __RESERVOIR_H__
#define __RESERVOIR_H__

#include "matrix.h"
#include "grid.h"

// Class for reservoir problem
class ReservoirProblem
{
   public:
      ReservoirProblem () {};
      ~ReservoirProblem () {};
      void run ();

   private:
      unsigned int max_iter;
      unsigned int nrk;
      double  ark[3], brk[3];
      double  cfl, final_time, dt;
      double  max_velocity;
      Grid    grid;
      Matrix  saturation;
      Matrix  concentration;
      Matrix  pressure;
      Matrix  permeability;

      void make_grid ();
      void initialize ();
      void residual (Matrix&, Matrix&);
      void solve ();
      void output (const unsigned int) const;

      std::vector<double> reconstruct
       (
       const unsigned int ill,
       const unsigned int jll,
       const unsigned int il,
       const unsigned int jl,
       const unsigned int ir,
       const unsigned int jr
       ) const;

      double darcy_velocity
         (const unsigned int, const unsigned int,
          const unsigned int, const unsigned int);

      void updateConcentration (Matrix&);
      void updateGhostCells ();
      void findMinMax () const;

};

double minmod (const double ul, const double u0, const double ur);

std::vector<double> num_flux
       (
       const double velocity,
       const std::vector<double> state_left,
       const std::vector<double> state_right
       );

#endif
