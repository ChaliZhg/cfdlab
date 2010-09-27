#include "matrix.h"
#include "grid.h"

class PressureProblem
{
   public:
      PressureProblem () {};
      PressureProblem (Grid*);
      ~PressureProblem () {};
      void run (const Matrix& saturation, 
                const Matrix& concentration,
                      Matrix& pressure);

   private:
      Grid*  grid;

      void compute_rhs ();
      void A_times_pressure (const Matrix& saturation,
                             const Matrix& concentration,
                             const Matrix& pressure);
};

