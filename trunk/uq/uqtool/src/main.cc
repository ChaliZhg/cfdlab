#include <iostream>
#include <cstdlib>
#include "uq.h"

using namespace std;

int main ()
{
   // Delete old RESULT directory
   system("rm -rf RESULT && mkdir RESULT");
   
   // Define UQ problem
   UQProblem<1> uq_problem;

   // Run UQ computation
   uq_problem.run ();
}
