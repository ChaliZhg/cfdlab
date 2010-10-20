#ifndef __PRESSURE_H__
#define __PRESSURE_H__

#include <vector>
#include <valarray>
#include "matrix.h"
#include "grid.h"

const double pinlet  = 1.0;
const double poutlet = 0.0;

class PressureProblem
{
   public:
      PressureProblem () {};
      PressureProblem (Grid*);
      ~PressureProblem () {};
      void run (const Matrix& saturation, 
                const Matrix& concentration,
                const Matrix& permeability,
                      Matrix& pressure);

   private:
      Grid*  grid;

      Matrix compute_rhs (const Matrix& saturation,
                          const Matrix& concentration,
                          const Matrix& permeability,
                          const Matrix& pressure);
      Matrix A_times_pressure (const Matrix& saturation,
                               const Matrix& concentration,
                               const Matrix& permeability,
                               const Matrix& pressure);
      Matrix residual (const Matrix& saturation,
                       const Matrix& concentration,
                       const Matrix& permeability,
                       const Matrix& pressure);
};

#endif
