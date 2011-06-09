#include <iostream>
#include <cstdlib>
#include "uq.h"

using namespace std;

// Run 1-variable uq
void problem_1d ()
{
   // Define UQ problem
   UQProblem<1> uq_problem;
   
   // Run UQ computation
   uq_problem.run ();
}

// Run 2-variable uq
void problem_2d ()
{
   // Define UQ problem
   UQProblem<2> uq_problem;
   
   // Run UQ computation
   uq_problem.run ();
}

// Main function
int main (int argc, char* argv[])
{
   if(argc != 2)
   {
      cout << "uqtool: specify problem dimension\n";
      abort ();
   }
   
   // Delete old RESULT directory
   system("rm -rf RESULT && mkdir RESULT");
   
   if(strcmp(argv[1],"1")==0)
      problem_1d ();
   else if(strcmp(argv[1],"2")==0)
      problem_2d ();
   else
   {
      cout << "Does not work for dimension = " << argv[1] << endl;
      abort ();
   }

  
}
