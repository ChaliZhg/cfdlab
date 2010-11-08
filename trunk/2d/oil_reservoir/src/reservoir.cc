#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include "pressure.h"
#include "reservoir.h"
#include "material.h"

#define SIGN(a) (((a)<0) ? -1:1)

using namespace std;

double minmod (const double ul, const double u0, const double ur)
{
   double result;

   double db = u0 - ul;         // backward difference
   double df = ur - u0;         // forward difference
   double dc = 0.5 * (ur - ul); // central difference

   if (db*dc > 0.0 && dc*df > 0.0)
   {
      result = min( min(fabs(db), fabs(dc)), fabs(df) );
      result *= SIGN(db);
   }
   else
      result = 0.0;

   return result;
}

// Find left state at interface between (il,jl) and (ir,jr)
// (ill,jll) is to the left of (il,jl)
vector<double> ReservoirProblem::reconstruct
       (
       const unsigned int ill,
       const unsigned int jll,
       const unsigned int il,
       const unsigned int jl,
       const unsigned int ir,
       const unsigned int jr
       ) const
{
   vector<double> state(2);
   double ds;

   // saturation
   ds = minmod (saturation(ill,jll), 
                saturation(il,jl), 
                saturation(ir,jr));
   state[0] = saturation (il,jl) + 0.5 * ds;

   // concentration
   // We reconstruct b = s*c
   double bll = saturation(ill,jll) * concentration(ill,jll);
   double bl  = saturation(il ,jl ) * concentration(il ,jl );
   double br  = saturation(ir ,jr ) * concentration(ir ,jr );
   double db  = minmod (bll, bl, br);
   state[1]   = (bl + 0.5 * db);
   if (state[0] > SZERO)
      state[1] /= state[0];
   else
      state[1] = 0.0;

   return state;
}

// computes total darcy velocity at interface b/w
// (ileft,jleft) and (iright,jright)
double ReservoirProblem::darcy_velocity
       (
       const unsigned int ileft,
       const unsigned int jleft,
       const unsigned int iright,
       const unsigned int jright
       )
{
   double s_left  = saturation    (ileft,  jleft);
   double c_left  = concentration (ileft,  jleft);
   double p_left  = pressure      (ileft,  jleft);
   double s_right = saturation    (iright, jright);
   double c_right = concentration (iright, jright);
   double p_right = pressure      (iright, jright);

   double m_total_left = mobility_total (s_left, c_left);
   double perm_left    = permeability (ileft, jleft);

   double m_total_right = mobility_total (s_right, c_right);
   double perm_right   = permeability (iright, jright);

   double m_perm   = harmonic_average (m_total_left  * perm_left, 
                                       m_total_right * perm_right);

   // we assume dx=dy, so dont divide by length since it is accounted
   // for in flux computation: dont multiply flux by face area
   double dpdn     = (p_right - p_left);
   double velocity = - m_perm * dpdn;

   max_velocity = max ( max_velocity, fabs(velocity)/grid.dx );
   min_velocity = min ( min_velocity, fabs(velocity)/grid.dx );
   
   return velocity;
}

// numerical flux function for saturation
vector<double> num_flux
       (
       const double velocity,
       const vector<double> state_left,
       const vector<double> state_right
       )
{
   double s_left  = state_left[0];
   double c_left  = state_left[1];
   double s_right = state_right[0];
   double c_right = state_right[1];

   double m_water_left = mobility_water (s_left, c_left);
   double m_oil_left   = mobility_oil (s_left, c_left);
   double m_total_left = m_water_left + m_oil_left;

   double m_water_right = mobility_water (s_right, c_right);
   double m_oil_right   = mobility_oil (s_right, c_right);
   double m_total_right = m_water_right + m_oil_right;

   vector<double> flux(2);

   if (velocity > 0)
   {
      flux[0] = velocity * m_water_left / m_total_left;
      flux[1] = c_left * flux[0];
   }
   else
   {
      flux[0] = velocity * m_water_right / m_total_right;
      flux[1] = c_right * flux[0];
   }

   return flux;
}

// Read some input from file
void ReservoirProblem::read_input ()
{
   cout << "Reading input from file data.in ..." << endl;

   ifstream inp;
   inp.open ("data.in");

   inp >> max_iter;
   inp >> cfl;
   inp >> cinlet;

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

   cout << "Max no. of time steps  = " << max_iter << endl;
   cout << "CFL number             = " << cfl << endl;
   cout << "nx x ny                = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries   = " << grid.n_boundary << endl;
   cout << "Number of cells        = " << grid.n_cells << endl;
   cout << "Number of actual cells = " << (grid.nx-1)*(grid.ny-1) << endl;

}

void ReservoirProblem::make_grid ()
{
   cout << "Making grid for reservoir problem ..." << endl;

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

   // This is implicitly assumed in the flux computations
   assert (grid.dx == grid.dy);

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

}

// allocate memory and set initial condition
void ReservoirProblem::initialize ()
{
   saturation.allocate    (grid.nx+1, grid.ny+1);
   concentration.allocate (grid.nx+1, grid.ny+1);
   pressure.allocate      (grid.nx+1, grid.ny+1);
   permeability.allocate  (grid.nx+1, grid.ny+1);

   // rock permeability
   for(unsigned int i=0; i<grid.nx+1; ++i)
      for(unsigned int j=0; j<grid.ny+1; ++j)
         permeability (i,j) = rock_permeability (grid.xc(i,j), grid.yc(i,j));

   // initialize only real cells, not for ghost cells
   for(unsigned int i=1; i<=grid.nx-1; ++i)
      for(unsigned int j=1; j<=grid.ny-1; ++j)
      {
         double dist = grid.xc(i,j) * grid.xc(i,j) + 
                       grid.yc(i,j) * grid.yc(i,j);
         if(dist <= 0.25 * 0.25)
         {
            saturation    (i,j) = 1.0;
            concentration (i,j) = cinlet;
         }
         else
         {
            saturation    (i,j) = 0.0;
            concentration (i,j) = 0.0;
         }

         // pressure is everywhere zero to begin with
         pressure (i,j) = 0.0;
      }

   updateGhostCells ();
   output (0);

   final_time = 0.1;
   nrk        = 3;
   ark[0] = 0.0; ark[1] = 3.0/4.0; ark[2] = 1.0/3.0;
   for (unsigned int i=0; i<3; ++i) brk[i] = 1.0 - ark[i];
}

// residual for saturation/concentration equation
void ReservoirProblem::residual (Matrix& s_residual, Matrix& c_residual)
{
   unsigned int i, j;
   double velocity;
   vector<double> state_left(2), state_right(2), flux(2);

   s_residual = 0.0;
   c_residual = 0.0;

   min_velocity = 1.0e20;
   max_velocity = 0.0;

   // interior vertical faces
   for(i=2; i<=grid.nx-1; ++i)
      for(j=1; j<=grid.ny-1; ++j)
      {
         state_left  = reconstruct (i-2, j, i-1, j, i, j);
         state_right = reconstruct (i+1, j, i, j, i-1, j);
         velocity    = darcy_velocity (i-1, j, i, j);
         flux        = num_flux (velocity, state_left, state_right);

         s_residual (i-1,j) += flux[0];
         s_residual (i,  j) -= flux[0];

         c_residual (i-1,j) += flux[1];
         c_residual (i,  j) -= flux[1];
      }

   // interior horizontal faces
   for(j=2; j<=grid.ny-1; ++j)
      for(i=1; i<=grid.nx-1; ++i)
      {
         state_left  = reconstruct (i, j+1, i, j, i, j-1);
         state_right = reconstruct (i, j-2, i, j-1, i, j);
         velocity    = darcy_velocity (i, j, i, j-1);
         flux        = num_flux (velocity, state_left, state_right);

         s_residual (i,j)   += flux[0];
         s_residual (i,j-1) -= flux[0];

         c_residual (i,j)   += flux[1];
         c_residual (i,j-1) -= flux[1];
      }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if (grid.ibeg[n] == grid.iend[n])
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {

            if (grid.ibeg[n] == 1) // inlet-vertical side
            {
               state_left  = reconstruct (i-1, j, i-1, j, i, j);
               state_right = reconstruct (i+1, j, i, j, i-1, j);
               velocity    = darcy_velocity (i-1, j, i, j);
               flux        = num_flux (velocity, state_left, state_right);
               s_residual(i,j) -= flux[0];
               c_residual(i,j) -= flux[1];
            }
            else // outlet-vertical side
            {
               state_left  = reconstruct (i-2, j, i-1, j, i, j);
               state_right = reconstruct (i, j, i, j, i-1, j);
               velocity    = darcy_velocity (i-1, j, i, j);
               flux        = num_flux (velocity, state_left, state_right);
               s_residual(i-1,j) += flux[0];
               c_residual(i-1,j) += flux[1];
            }
         }
      }

      if (grid.jbeg[n] == grid.jend[n])
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {

            if(grid.jbeg[n] == 1) // inlet-horizontal side
            {
               state_left  = reconstruct (i, j+1, i, j, i, j-1);
               state_right = reconstruct (i, j-1, i, j-1, i, j);
               velocity    = darcy_velocity (i, j, i, j-1);
               flux        = num_flux (velocity, state_left, state_right);
               s_residual(i,j) += flux[0];
               c_residual(i,j) += flux[1];
            }
            else // outlet-horizontal side
            {
               state_left  = reconstruct (i, j, i, j, i, j-1);
               state_right = reconstruct (i, j-2, i, j-1, i, j);
               velocity    = darcy_velocity (i, j, i, j-1);
               flux        = num_flux (velocity, state_left, state_right);
               s_residual(i,j-1) -= flux[0];
               c_residual(i,j-1) -= flux[1];
            }
         }
      }
   }

   dt = cfl * max (grid.dx, grid.dy) / (3.0 * max_velocity);
   double lambda = dt / (grid.dx * grid.dy);
   s_residual *= lambda;
   c_residual *= lambda;
}

// Update polymer concentration
void ReservoirProblem::updateConcentration (Matrix& sc)
{
   for (unsigned int i=1; i<grid.nx; ++i)
      for (unsigned int j=1; j<grid.ny; ++j)
      {
         if (saturation (i,j) > SZERO)
            concentration (i,j) = sc (i,j) / saturation (i,j);
         else
            concentration (i,j) = 0.0;
      }
}

// Update solution in ghost cells
void ReservoirProblem::updateGhostCells ()
{
   unsigned int i, j;

   // top/bottom ghost cells
   for (i=1; i<grid.nx; ++i)
   {
      j = 0;
      saturation    (i,j) = saturation    (i,j+1);
      concentration (i,j) = concentration (i,j+1);
      pressure      (i,j) = pressure      (i,j+1);

      j = grid.ny;
      saturation    (i,j) = saturation    (i,j-1);
      concentration (i,j) = concentration (i,j-1);
      pressure      (i,j) = pressure      (i,j-1);
   }

   // left/right ghost cells
   for (j=1; j<grid.ny; ++j)
   {
      i = 0;
      saturation    (i,j) = saturation    (i+1,j);
      concentration (i,j) = concentration (i+1,j);
      pressure      (i,j) = pressure      (i+1,j);

      i = grid.nx;
      saturation    (i,j) = saturation    (i-1,j);
      concentration (i,j) = concentration (i-1,j);
      pressure      (i,j) = pressure      (i-1,j);
   }

   // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if (grid.ibeg[n] == grid.iend[n])
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {

            if (grid.ibeg[n] == 1) // inlet-vertical side
            {
               saturation    (i-1,j) = 1.0;
               concentration (i-1,j) = cinlet;
               pressure      (i-1,j) = pinlet;
            }
            else // outlet-vertical side
            {
               //saturation    (i,j) = 0.0;
               //concentration (i,j) = 0.0;
               pressure      (i,j) = poutlet;
            }
         }
      }

      if (grid.jbeg[n] == grid.jend[n])
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {

            if(grid.jbeg[n] == 1) // inlet-horizontal side
            {
               saturation    (i,j-1) = 1.0;
               concentration (i,j-1) = cinlet;
               pressure      (i,j-1) = pinlet;
            }
            else // outlet-horizontal side
            {
               //saturation    (i,j) = 1.0;
               //concentration (i,j) = 0.0;
               pressure      (i,j) = poutlet;
            }
         }
      }
   }
}

// Find min and max values of solution
void ReservoirProblem::findMinMax () const
{
   double s_min = 1.0e20;
   double s_max =-1.0e20;
   double c_min = 1.0e20;
   double c_max =-1.0e20;
   double p_min = 1.0e20;
   double p_max =-1.0e20;

   for (unsigned int i=1; i<grid.nx; ++i)
      for (unsigned int j=1; j<grid.ny; ++j)
      {
         s_min = min (s_min, saturation(i,j));
         s_max = max (s_max, saturation(i,j));
         c_min = min (c_min, concentration(i,j));
         c_max = max (c_max, concentration(i,j));
         p_min = min (p_min, pressure(i,j));
         p_max = max (p_max, pressure(i,j));
      }

   cout << "Saturation    = " << s_min << " " << s_max << endl;
   cout << "Concentration = " << c_min << " " << c_max << endl;
   cout << "Pressure      = " << p_min << " " << p_max << endl;
   cout << "Min velocity  = " << min_velocity << endl;
   cout << "Max velocity  = " << max_velocity << endl;
   cout << "dt            = " << dt << endl;
}

// perform time stepping
void ReservoirProblem::solve ()
{
   unsigned int iter = 0;
   double time = 0.0;
   PressureProblem pressure_problem (&grid);
   Matrix s_residual (grid.nx+1, grid.ny+1);
   Matrix c_residual (grid.nx+1, grid.ny+1);
   Matrix sc         (grid.nx+1, grid.ny+1);
   Matrix s_old      (grid.nx+1, grid.ny+1);
   Matrix sc_old     (grid.nx+1, grid.ny+1);

   while (iter < max_iter)
   {
      s_old = saturation;
      sc_old= saturation * concentration;

      // solve for pressure
      pressure_problem.run (saturation, concentration, 
                            permeability, pressure);

      for (unsigned int irk=0; irk<nrk; ++irk)
      {
         // compute residual
         residual (s_residual, c_residual);
      
         // new value of s*c
         sc = sc_old * ark[irk] + 
              (saturation * concentration - c_residual) * brk[irk];

         // update saturation
         saturation  = s_old * ark[irk] + (saturation - s_residual) * brk[irk];

         // update concentration
         updateConcentration (sc);

         // update solution in ghost cells
         updateGhostCells ();
      }

      // find solution range: to check for stability
      findMinMax ();
      
      time += dt;
      ++iter;

      // save solution to file
      if (iter % 25 == 0)
         output (iter);

      cout << "Time= " << time << " iter= " << iter << endl;
      cout << endl;
   }
}

// save solution to file
// only interior cells are written, ghost cells are not written
void ReservoirProblem::output (const unsigned int iter) const
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
         vtk << fixed << pressure (i,j) << endl;

   // write saturation
   vtk << "SCALARS saturation float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << saturation (i,j) << endl;

   // write concentration
   vtk << "SCALARS concentration float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << concentration (i,j) << endl;

   // write permeability
   if (iter==0)
   {
      vtk << "SCALARS permeability float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(j=1; j<grid.ny; ++j)
         for(i=1; i<grid.nx; ++i)
            vtk << fixed << permeability (i,j) << endl;
   }

   vtk.close ();

}

// solve the whole problem
void ReservoirProblem::run ()
{
   read_input ();
   make_grid ();
   initialize ();
   solve ();
}
