#include <iostream>
#include <cstdlib>
#include "uq.h"

using namespace std;

int main ()
{
   // Delete old sample directories S000-S999
   system("rm -rf RESULT && mkdir RESULT");
   
   // Define UQ problem
   UQProblem<1> uq_problem;

   // Run UQ computation
   uq_problem.run ();
}
