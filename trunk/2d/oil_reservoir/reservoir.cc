#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
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

   for(unsigned int n=0; n<grid.n_boundary; ++n){
      inp >> grid.ibeg[n] >> grid.iend[n]
          >> grid.jbeg[n] >> grid.jend[n]
          >> grid.boundary_condition[n];
      assert (grid.ibeg[n] >= 1 && grid.ibeg[n] <= grid.nx);
      assert (grid.iend[n] >= 1 && grid.iend[n] <= grid.nx);
      assert (grid.jbeg[n] >= 1 && grid.jbeg[n] <= grid.ny);
      assert (grid.jend[n] >= 1 && grid.jend[n] <= grid.ny);
      if(grid.ibeg[n] == grid.iend[n] &&
         grid.jbeg[n] == grid.jend[n])
      {
         cout << "Boundary " << n 
              << " is not a surface !!!" << endl;
         abort ();
      }
   }

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

   for(unsigned int i=1; i<=grid.nx; ++i)
      for(unsigned int j=1; j<=grid.ny; ++j)
      {
         grid.x (i,j) = xmin + (i-1) * grid.dx;
         grid.y (i,j) = ymin + (j-1) * grid.dy;
      }

   cout << "nx x ny                = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries   = " << grid.n_boundary << endl;
   cout << "Number of cells        = " << grid.n_cells << endl;
   cout << "Number of actual cells = " << (grid.nx-1)*(grid.ny-1) << endl;

}

// allocate memory and initial condition
void ReservoirProblem::initialize ()
{
   saturation.allocate    (grid.nx+1, grid.ny+1);
   concentration.allocate (grid.nx+1, grid.ny+1);
   pressure.allocate      (grid.nx+1, grid.ny+1);

   for(unsigned int i=1; i<=grid.nx-1; ++i)
      for(unsigned int j=1; j<=grid.ny-1; ++j)
      {
         saturation    (i,j) = 0.0;
         concentration (i,j) = 0.0;
         pressure      (i,j) = 0.0;
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

      if (iter == 1 || iter % 10 == 0)
         output (iter);

      cout << "Time= " << time << " iter= " << iter << endl;
   }
}

// save solution to file
void ReservoirProblem::output (const unsigned int iter)
{

   unsigned int i, j;
   ofstream vtk;
   ostringstream filename;
   filename << "solution-" << iter << ".vtk";

   vtk.open (filename.str().c_str());

   vtk << "# vtk DataFile Version 2.0" << endl;
   vtk << "Oil Reservoir Problem: iter = " << iter << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_GRID" << endl;
   vtk << "DIMENSIONS " << grid.nx << " " << grid.ny << " 1" << endl;

   // write coordinates
   vtk << "POINTS " << grid.nx * grid.ny << " float" << endl;
   for(j=1; j<=grid.ny; ++j)
      for(i=1; i<=grid.nx; ++i)
      {
         vtk << grid.x (i,j) << "  "
             << grid.y (i,j) << "  "
             << 0.0 << endl;
      }

   vtk << "CELL_DATA " << (grid.nx-1)*(grid.ny-1) << endl;

   // write pressure
   vtk << "SCALARS pressure float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << pressure (i,j) << endl;

   vtk.close ();

}

// solve the whole problem
void ReservoirProblem::run ()
{
   make_grid ();
   initialize ();
   solve ();
}
