#include <iostream>
#include <cstdlib>
#include "uq.h"

using namespace std;

int main ()
{
   // Delete old sample directories S000-S999
   system("rm -rf S[0-9][0-9][0-9]");
   system("rm -rf EVAL");
   
   // Define UQ problem
   UQProblem<1> uq_problem;

   // Run UQ computation
   uq_problem.run ();
}
