#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include "pressure.h"
#include "reservoir.h"
#include "material.h"

using namespace std;

// numerical flux function for saturation
double s_num_flux ()
{
   return 1.0;
}

// numerical flux function for concentration
double c_num_flux ()
{
   return 1.0;
}

// Read some input and make the grid
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

   // grid vertex coordinates
   for(unsigned int i=0; i<=grid.nx+1; ++i)
      for(unsigned int j=0; j<=grid.ny+1; ++j)
      {
         grid.x (i,j) = xmin + (i-1) * grid.dx;
         grid.y (i,j) = ymin + (j-1) * grid.dy;
      }

   // cell center coordinates
   for(unsigned int i=0; i<=grid.nx; ++i)
      for(unsigned int j=0; j<=grid.ny; ++j)
      {
         grid.xc (i,j) = 0.25 * ( grid.x(i,j)     + grid.x(i+1,j) + 
                                  grid.x(i+1,j+1) + grid.x(i,j+1) );
         grid.yc (i,j) = 0.25 * ( grid.y(i,j)     + grid.y(i+1,j) + 
                                  grid.y(i+1,j+1) + grid.y(i,j+1) );
      }

   cout << "nx x ny                = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries   = " << grid.n_boundary << endl;
   cout << "Number of cells        = " << grid.n_cells << endl;
   cout << "Number of actual cells = " << (grid.nx-1)*(grid.ny-1) << endl;

}

// allocate memory and set initial condition
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

// residual for saturation/concentration equation
void ReservoirProblem::residual (Matrix& s_residual, Matrix& c_residual)
{
   unsigned int i, j;
   double mobility_left, mobility_right, mobility;
   double dpdn, velocity, flux;

   s_residual = 0.0;
   c_residual = 0.0;

   // interior vertical faces
   for(i=2; i<=grid.nx-1; ++i)
      for(j=1; j<=grid.ny-1; ++j)
      {
         mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
         mobility_right = mobility_total (saturation(i,j), concentration(i,j));
         mobility = harmonic_average (mobility_left, mobility_right);

         dpdn     = (pressure(i,j) - pressure(i-1,j))/grid.dx;
         velocity = - mobility * dpdn;

         flux                = s_num_flux () * grid.dy;
         s_residual (i-1,j) += flux;
         s_residual (i,  j) -= flux;

         flux                = c_num_flux () * grid.dy;
         c_residual (i-1,j) += flux;
         c_residual (i,  j) -= flux;
      }

   // interior horizontal faces
   for(j=2; j<=grid.ny-1; ++j)
      for(i=1; i<=grid.nx-1; ++i)
      {
         mobility_left = mobility_total (saturation(i,j), concentration(i,j));
         mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
         mobility = harmonic_average (mobility_left, mobility_right);

         dpdn     = (pressure(i,j-1) - pressure(i,j))/grid.dy;
         velocity = - mobility * dpdn;

         flux                = s_num_flux () * grid.dx;
         s_residual (i,j)   += flux;
         s_residual (i,j-1) -= flux;

         flux                = c_num_flux () * grid.dx;
         c_residual (i,j)   += flux;
         c_residual (i,j-1) -= flux;
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if (grid.ibeg[n] == grid.iend[n])
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {
            mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            mobility = harmonic_average (mobility_left, mobility_right);

            if (grid.ibeg[n] == 1) // inlet-vertical side
            {
               dpdn     = (pressure(i,j) - pinlet)/(0.5 * grid.dx);
               velocity = - mobility * dpdn;

               flux             = s_num_flux () * grid.dy;
               s_residual(i,j) -= flux;

               flux             = c_num_flux () * grid.dy;
               c_residual(i,j) -= flux;
            }
            else // outlet-vertical side
            {
               dpdn     = (poutlet - pressure(i-1,j))/(0.5 * grid.dx);
               velocity = - mobility * dpdn;

               flux               = s_num_flux () * grid.dy;
               s_residual(i-1,j) += flux;

               flux               = c_num_flux () * grid.dy;
               c_residual(i-1,j) += flux;
            }
         }
      }

      if (grid.jbeg[n] == grid.jend[n])
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {
            mobility_left = mobility_total (saturation(i,j), concentration(i,j));
            mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
            mobility = harmonic_average (mobility_left, mobility_right);

            if(grid.jbeg[n] == 1) // inlet-horizontal side
            {
               dpdn     = (pinlet - pressure(i,j))/(0.5 * grid.dy);
               velocity = - mobility * dpdn;

               flux             = s_num_flux () * grid.dx;
               s_residual(i,j) += flux;

               flux             = c_num_flux () * grid.dx;
               c_residual(i,j) += flux;
            }
            else // outlet-horizontal side
            {
               dpdn     = (pressure(i,j-1) - poutlet)/(0.5 * grid.dy);
               velocity = - mobility * dpdn;

               flux               = s_num_flux () * grid.dx;
               s_residual(i,j-1) -= flux;

               flux               = c_num_flux () * grid.dx;
               c_residual(i,j-1) -= flux;
            }
         }
      }
   }

   double lambda = dt / (grid.dx * grid.dy);
   s_residual *= lambda;
   c_residual *= lambda;
}

// perform time stepping
void ReservoirProblem::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   PressureProblem pressure_problem (&grid);
   Matrix s_residual (grid.nx+1, grid.ny+1);
   Matrix c_residual (grid.nx+1, grid.ny+1);

   while (time < final_time)
   {
      // solve for pressure
      pressure_problem.run (saturation, concentration, pressure);

      // compute residual
      residual (s_residual, c_residual);
      
      // update saturation/concentration
      saturation    -= s_residual;
      concentration -= c_residual;

      // update solution in ghost cells
      
      time += dt;
      ++iter;

      // save solution to file
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

   // write saturation
   vtk << "SCALARS saturation float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << saturation (i,j) << endl;

   // write concentration
   vtk << "SCALARS concentration float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << concentration (i,j) << endl;

   vtk.close ();

}

// solve the whole problem
void ReservoirProblem::run ()
{
   make_grid ();
   initialize ();
   solve ();
}
