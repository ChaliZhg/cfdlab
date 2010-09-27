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
      double  final_time, dt;
      Grid    grid;
      Matrix  saturation;
      Matrix  concentration;
      Matrix  pressure;

      void make_grid ();
      void initialize ();
      void solve ();
      void output ();

};
