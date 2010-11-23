#include <iostream>
#include "reservoir.h"

double gravity;

using namespace std;

int main ()
{
   cout << "Starting reservoir problem ..." << endl;

   ReservoirProblem  reservoir_problem;

   reservoir_problem.run ();

   return 0;
}
