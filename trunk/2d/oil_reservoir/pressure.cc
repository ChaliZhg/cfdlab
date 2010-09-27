#include <iostream>
#include "matrix.h"
#include "grid.h"
#include "pressure.h"
#include "material.h"

#define harmonic_average(a,b)   (2.0*(a)*(b)/((a)+(b)))

using namespace std;

// constructor given grid
PressureProblem::PressureProblem (Grid* grid_in)
{
   grid = grid_in;
}

// compute rhs (b) of pressure equation A*p=b
void PressureProblem::compute_rhs ()
{
}

// compute matrix vector product A*b in pressure equation A*p = b
void PressureProblem::A_times_pressure (const Matrix& saturation,
                                        const Matrix& concentration,
                                        const Matrix& pressure)
{
   unsigned int i, j, left_cell, right_cell;
   double mobility, mobility_left, mobility_right;
   double dpdn, flux;
   vector<double> result(grid->n_cells, 0.0);

   // interior vertical faces
   for(i=1; i<grid->nx-1; ++i)
      for(j=0; j<grid->ny-1; ++j)
      {
         mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
         mobility_right = mobility_total (saturation(i,j), concentration(i,j));
         mobility = harmonic_average (mobility_left, mobility_right);

         dpdn = (pressure(i,j) - pressure(i-1,j))/grid->dx;

         flux                = mobility * dpdn;

         left_cell           = grid->cell_num (i-1, j);
         result[left_cell]  += flux;

         right_cell          = grid->cell_num (i,j);
         result[right_cell] -= flux;
      }

   // interior horizontal faces
   for(j=1; j<grid->ny-1; ++j)
      for(i=0; i<grid->nx-1; ++i)
      {
         mobility_left = mobility_total (saturation(i,j), concentration(i,j));
         mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
         mobility = harmonic_average (mobility_left, mobility_right);

         dpdn = (pressure(i,j-1) - pressure(i,j))/grid->dy;

         flux                = mobility * dpdn;

         left_cell           = grid->cell_num (i, j);
         result[left_cell]  += flux;

         right_cell          = grid->cell_num (i,j-1);
         result[right_cell] -= flux;
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid->n_boundary; ++n)
   {
      for(i=grid->ibeg[n]; i<grid->iend[n]; ++i)
         for(j=grid->jbeg[n]; j<grid->jend[n]; ++j)
         {
            if(grid->b_type[n] == imin)
            {
            }
         }
   }

}

// Solve pressure equation by CG method
void PressureProblem::run (const Matrix& saturation, 
                           const Matrix& concentration,
                                 Matrix& pressure)
{
   const unsigned int max_iter = 1000;
   const double tolerance = 1.0e-6;
   unsigned int iter = 0;
   double residue = 1.0;
   vector<double> direction(grid->n_cells);

   compute_rhs ();

   while ( residue > tolerance && iter < max_iter )
   {
      A_times_pressure (saturation, concentration, pressure);
      ++iter;
   }

   if (residue > tolerance && iter==max_iter)
   {
      cout << "PressureProblem did not converge !!!" << endl;
   }
}
