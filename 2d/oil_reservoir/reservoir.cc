#include <iostream>
#include <fstream>
#include "pressure.h"
#include "reservoir.h"

using namespace std;

void ReservoirProblem::make_grid ()
{
   cout << "Making grid for reservoir problem ..." << endl;

   ifstream inp;
   inp.open ("data.in");

   inp >> grid.nx >> grid.ny;
   inp >> grid.n_boundary;

   // allocate memory for grid
   grid.allocate ();

   for(unsigned int n=0; n<grid.n_boundary; ++n)
      inp >> grid.ibeg[n] >> grid.iend[n]
          >> grid.jbeg[n] >> grid.jend[n]
          >> grid.boundary_condition[n];

   inp.close ();

   // set location of boundary
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if(grid.ibeg[n] == grid.iend[n])
      {
         if(grid.ibeg[n] == 0)
            grid.b_type[n] = imin;
         else
            grid.b_type[n] = imax;
      }

      if(grid.jbeg[n] == grid.jend[n])
      {
         if(grid.jbeg[n] == 0)
            grid.b_type[n] = jmin;
         else
            grid.b_type[n] = jmax;
      }
   }

   double xmin = 0.0;
   double xmax = 1.0;
   double ymin = 0.0;
   double ymax = 1.0;
   grid.dx   = (xmax - xmin)/(grid.nx - 1);
   grid.dy   = (ymax - ymin)/(grid.ny - 1);

   for(unsigned int i=0; i<grid.nx; ++i)
      for(unsigned int j=0; j<grid.ny; ++j)
      {
         grid.x (i,j) = xmin + i*grid.dx;
         grid.y (i,j) = ymin + j*grid.dy;
      }

   cout << "nx x ny              = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries = " << grid.n_boundary << endl;
   cout << "Number of cells      = " << grid.n_cells << endl;

}

// allocate memory and initial condition
void ReservoirProblem::initialize ()
{
   saturation.allocate (grid.nx-1, grid.ny-1);
   concentration.allocate (grid.nx-1, grid.ny-1);
   pressure.allocate (grid.nx-1, grid.ny-1);

   for(unsigned int i=0; i<grid.nx-1; ++i)
      for(unsigned int j=0; j<grid.ny-1; ++j)
      {
         saturation(i,j) = 0.0;
         concentration(i,j) = 0.0;
         pressure(i,j) = 0.0;
      }

   final_time = 0.01;
   dt         = 0.01;
}

// perform time stepping
void ReservoirProblem::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   PressureProblem pressure_problem (&grid);

   while (time < final_time)
   {
      // solve for pressure
      pressure_problem.run (saturation, concentration, pressure);

      // compute residual
      
      // update saturation/concentration
      
      time += dt;
      ++iter;

      cout << "Time= " << time << " iter= " << iter << endl;
   }
}

// save solution to file
void ReservoirProblem::output ()
{

   ofstream vtk;
   vtk.open ("solution.vtk");
   for(unsigned int i=0; i<grid.nx; ++i)
      for(unsigned int j=0; j<grid.ny; ++j)
      {
         vtk << grid.x (i,j) << "  ";
         vtk << grid.y (i,j) << endl;
      }
   vtk.close ();

}

// solve the whole problem
void ReservoirProblem::run ()
{
   make_grid ();
   initialize ();
   solve ();
   output ();
}
