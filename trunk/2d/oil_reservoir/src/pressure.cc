#include <iostream>
#include <valarray>
#include "matrix.h"
#include "grid.h"
#include "pressure.h"
#include "material.h"

using namespace std;

// constructor given grid
PressureProblem::PressureProblem (Grid* grid_in)
{
   // set pointer grid to grid_in
   grid = grid_in;
}

// compute rhs (b) of pressure equation A*p=b
Matrix PressureProblem::compute_rhs (const Matrix& saturation,
                                     const Matrix& concentration,
                                     const Matrix& permeability,
                                     const Matrix& pressure)
{
   unsigned int i, j;
   double mobility_left, mobility_right;
   double perm_left, perm_right;
   double m_perm;
   double flux;
   Matrix result(grid->nx+1, grid->ny+1);

   result = 0.0;

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid->n_boundary; ++n)
   {
      if (grid->ibeg[n] == grid->iend[n])
      {
         i = grid->ibeg[n];
         for(j=grid->jbeg[n]; j<grid->jend[n]; ++j)
         {
            mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
            perm_left = permeability (i-1,j);
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            perm_right = permeability (i,j);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if (grid->ibeg[n] == 1) // inlet-vertical side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dx)
               flux         = m_perm * (-pinlet)/(grid->dx) * grid->dy;
               result(i,j) -= flux;
            }
            else // outlet-vertical side
            {
               // dpdn = (poutlet - pressure(i-1,j))/(dx)
               flux           = m_perm * (poutlet)/(grid->dx) * grid->dy;
               result(i-1,j) += flux;
            }
         }
      }

      if (grid->jbeg[n] == grid->jend[n])
      {
         j = grid->jbeg[n];
         for(i=grid->ibeg[n]; i<grid->iend[n]; ++i)
         {
            mobility_left = mobility_total (saturation(i,j), concentration(i,j));
            perm_left = permeability (i,j);
            mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
            perm_right = permeability (i,j-1);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if(grid->jbeg[n] == 1) // inlet-horizontal side
            {
               // dpdn = (pinlet - pressure(i,j))/(dy)
               flux         = m_perm * (pinlet)/(grid->dy) * grid->dx;
               result(i,j) += flux;
            }
            else // outlet-horizontal side
            {
               // dpdn = (pressure(i,j-1) - poutlet)/(dy)
               flux           = m_perm * (-poutlet)/(grid->dy) * grid->dx;
               result(i,j-1) -= flux;
            }
         }
      }
   }

   return result;

}

// compute matrix vector product A*p in pressure equation A*p = b
Matrix PressureProblem::A_times_pressure (const Matrix& saturation,
                                          const Matrix& concentration,
                                          const Matrix& permeability,
                                          const Matrix& pressure)
{
   unsigned int i, j;
   double mobility_left, mobility_right;
   double perm_left, perm_right;
   double m_perm, dpdn, flux;
   Matrix result(grid->nx+1, grid->ny+1);

   result = 0.0;

   // interior vertical faces
   for(i=2; i<=grid->nx-1; ++i)
      for(j=1; j<=grid->ny-1; ++j)
      {
         mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
         perm_left = permeability (i-1,j);
         mobility_right = mobility_total (saturation(i,j), concentration(i,j));
         perm_right = permeability (i,j);
         m_perm = harmonic_average (mobility_left  * perm_left, 
                                    mobility_right * perm_right);

         dpdn = (pressure(i,j) - pressure(i-1,j))/grid->dx;

         flux           = m_perm * dpdn * grid->dy;
         result(i-1,j) += flux;
         result(i,j)   -= flux;
      }

   // interior horizontal faces
   for(j=2; j<=grid->ny-1; ++j)
      for(i=1; i<=grid->nx-1; ++i)
      {
         mobility_left = mobility_total (saturation(i,j), concentration(i,j));
         perm_left = permeability (i,j);
         mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
         perm_right = permeability (i,j-1);
         m_perm = harmonic_average (mobility_left  * perm_left, 
                                    mobility_right * perm_right);

         dpdn = (pressure(i,j-1) - pressure(i,j))/grid->dy;

         flux           = m_perm * dpdn * grid->dx;
         result(i,j)   += flux;
         result(i,j-1) -= flux;
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid->n_boundary; ++n)
   {
      if (grid->ibeg[n] == grid->iend[n])
      {
         i = grid->ibeg[n];
         for(j=grid->jbeg[n]; j<grid->jend[n]; ++j)
         {
            mobility_left = mobility_total (saturation(i-1,j), concentration(i-1,j));
            perm_left = permeability (i-1,j);
            mobility_right = mobility_total (saturation(i,j), concentration(i,j));
            perm_right = permeability (i,j);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if (grid->ibeg[n] == 1) // inlet-vertical side
            {
               // dpdn = (pressure(i,j) - pinlet)/(dx)
               flux         = m_perm * pressure(i,j)/(grid->dx) * grid->dy;
               result(i,j) -= flux;
            }
            else // outlet-vertical side
            {
               // dpdn = (poutlet - pressure(i-1,j))/(dx)
               flux           = m_perm * (-pressure(i-1,j))/(grid->dx) * grid->dy;
               result(i-1,j) += flux;
            }
         }
      }

      if (grid->jbeg[n] == grid->jend[n])
      {
         j = grid->jbeg[n];
         for(i=grid->ibeg[n]; i<grid->iend[n]; ++i)
         {
            mobility_left = mobility_total (saturation(i,j), concentration(i,j));
            perm_left = permeability (i,j);
            mobility_right = mobility_total (saturation(i,j-1), concentration(i,j-1));
            perm_right = permeability (i,j-1);
            m_perm = harmonic_average (mobility_left  * perm_left, 
                                       mobility_right * perm_right);

            if(grid->jbeg[n] == 1) // inlet-horizontal side
            {
               // dpdn = (pinlet - pressure(i,j))/(dy)
               flux         = m_perm * (-pressure(i,j))/(grid->dy) * grid->dx;
               result(i,j) += flux;
            }
            else // outlet-horizontal side
            {
               // dpdn = (pressure(i,j-1) - poutlet)/(dy)
               flux           = m_perm * pressure(i,j-1)/(grid->dy) * grid->dx;
               result(i,j-1) -= flux;
            }
         }
      }
   }

   // We need negative since -div(lambda*K*grad(p)) = source
   result *= -1.0;

   return result;

}

// Compute residual for pressure problem, r = b - A*p
Matrix PressureProblem::residual (const Matrix& saturation, 
                                  const Matrix& concentration,
                                  const Matrix& permeability,
                                  const Matrix& pressure)
{
   Matrix r(grid->nx+1, grid->ny+1);

   r = compute_rhs (saturation, concentration, permeability, pressure)
     - A_times_pressure(saturation, concentration, permeability, pressure);

   return r;
}

// Solve pressure equation by CG method
void PressureProblem::run (const Matrix& saturation, 
                           const Matrix& concentration,
                           const Matrix& permeability,
                                 Matrix& pressure)
{
   const unsigned int max_iter = 5000;
   const double tolerance = 1.0e-6;
   unsigned int iter = 0;
   double beta, omega;
   vector<double> r2(max_iter);
   Matrix d (grid->nx + 1, grid->ny + 1);
   Matrix r (grid->nx + 1, grid->ny + 1);
   Matrix v (grid->nx + 1, grid->ny + 1);

   // initial residual
   r = residual (saturation, concentration, permeability, pressure);

   // initial direction
   d = r;

   r2[0] = r.dot(r);

   // CG iterations
   while ( sqrt(r2[iter]) > tolerance && iter < max_iter )
   {
      if (iter >= 1) // update descent direction
      {              // d = r + beta * d
         beta = r2[iter] / r2[iter-1];
         d   *= beta;
         d   += r;
      }

      v = A_times_pressure (saturation, concentration, permeability, d);
      omega = r2[iter] / d.dot(v);

      // update pressure: p = p + omega * d
      pressure += d * omega;

      // update residual: r = r - omega * v
      v *= omega;
      r -= v;

      ++iter;

      r2[iter] = r.dot(r);

   }
   cout << "PressureProblem: iter= " << iter 
        << " residue= " << sqrt(r2[iter]) << endl;

   if (sqrt(r2[iter]) > tolerance && iter==max_iter)
   {
      cout << "PressureProblem did not converge !!!" << endl;
      abort ();
   }
}
