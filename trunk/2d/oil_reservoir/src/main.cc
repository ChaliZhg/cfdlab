#include <iostream>
#include "reservoir.h"

// Global variables
double gravity;
double pinlet;
double poutlet;
double viscosity_oil;
double density_water;
double density_oil;

using namespace std;

int main ()
{
   cout << "Starting reservoir problem ..." << endl;

   ReservoirProblem  reservoir_problem;

   reservoir_problem.run ();

   return 0;
}
