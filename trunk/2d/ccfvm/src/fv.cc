#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "parameter.h"
#include "fv.h"
#include "writer.h"

extern Dimension dim;
extern bool restart;
extern bool preprocess;
extern bool bounds;
extern bool convert_to_vtk;
extern bool convert_to_tec;

using namespace std;

bool check_for_stop_file ();

const double gauss_x[2] = {-1.0/sqrt(3.0), +1.0/sqrt(3.0)};
const double gauss_wt[2] = {0.5, 0.5};

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FiniteVolume::initialize ()
{
   cout << "Initializing memory\n";
   primitive.resize (grid.n_cell);
   conserved.resize (grid.n_cell);
   conserved_old.resize (grid.n_cell);
   ang_mom.resize (grid.n_cell);
   ang_mom_old.resize (grid.n_cell);
   residual.resize (grid.n_cell);
   res_ang_mom.resize (grid.n_cell);
   phi.resize (grid.n_cell);
   sdxdU.resize (grid.n_cell);
   sdydU.resize (grid.n_cell);
   if(param.smooth_res)
   {
      residual1.resize (grid.n_cell);
      residual2.resize (grid.n_cell);
   }
   dt.resize (grid.n_cell);

   // Gradients of conserved variables
   grad_rho.resize (grid.n_cell);
   grad_rhoU.resize (grid.n_cell);
   grad_rhoV.resize (grid.n_cell);
   grad_E.resize (grid.n_cell);

   // If restart option specified, read previous solution from file
   if(restart)
   {
      cout << "Reading restart file restart.dat ...\n";
      ifstream fi;
      fi.open ("restart.dat");
      assert (fi.is_open());
      for(unsigned int i=0; i<grid.n_cell; ++i)
         fi >> primitive[i].density
            >> primitive[i].velocity.x
            >> primitive[i].velocity.y
            >> primitive[i].pressure;
      fi >> last_iter;
      fi.close ();
   }
   else
   {
      cout << "Setting initial condition to input values ...";
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         primitive[i] = param.initial_condition.value (grid.cell[i].centroid);
         assert (primitive[i].density  > 0.0);
         assert (primitive[i].pressure > 0.0);
      }
      cout << " Done\n";
      last_iter = 0;
   }

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      conserved[i] = param.material.prim2con(primitive[i]);
      ang_mom[i]   = 0.0; // TODO
   }
   
   // Check if solution conversion was requested
   // Then save solution and abort
   if(convert_to_tec) param.write_format = "tec";
   if(convert_to_vtk) param.write_format = "vtk";
   if(convert_to_tec || convert_to_vtk)
   {
      cout << "Saving solution in ";
      if(convert_to_tec) cout << "tecplot ";
      if(convert_to_vtk) cout << "vtk ";
      cout << "format\n";
      compute_gradients ();
      output (0);
   }
}

//------------------------------------------------------------------------------
// Compute derivatives of conserved variables
//------------------------------------------------------------------------------
void FiniteVolume::compute_gradients ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      sdxdU[i] = 0.0;
      sdydU[i] = 0.0;
   }
   
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      int cl = grid.face[i].lcell;
      int cr = grid.face[i].rcell;
      
      ConVar dU;
      Vector dr;
      if(cr >= 0)
      {
         dU = conserved[cr] - conserved[cl];
         dr = grid.cell[cr].centroid - grid.cell[cl].centroid;
      }
      else // This is a boundary face
      {
         // apply boundary condition
         int face_type = grid.face[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         std::vector<ConVar> state(2);
         state[0] = conserved[cl];
         state[1] = conserved[cl];
         bc.apply (grid.face[i].centroid, grid.face[i], state);
         dU = state[1] - conserved[cl];
         dr = grid.face[i].centroid - grid.cell[cl].centroid;
      }
      
      sdxdU[cl] += dU * dr.x;
      sdydU[cl] += dU * dr.y;
      
      if(cr >= 0)
      {
         sdxdU[cr] += dU * dr.x;
         sdydU[cr] += dU * dr.y;
      }
   }
   
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      grad_rho[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].density +
                      grid.cell[i].invA1[0][1] * sdydU[i].density;
      grad_rho[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].density +
                      grid.cell[i].invA1[1][1] * sdydU[i].density;
      
      grad_rhoU[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].momentum.x +
                       grid.cell[i].invA1[0][1] * sdydU[i].momentum.x;
      grad_rhoU[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].momentum.x +
                       grid.cell[i].invA1[1][1] * sdydU[i].momentum.x;

      grad_rhoV[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].momentum.y +
                       grid.cell[i].invA1[0][1] * sdydU[i].momentum.y;
      grad_rhoV[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].momentum.y +
                       grid.cell[i].invA1[1][1] * sdydU[i].momentum.y;
      
      grad_E[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].energy +
                    grid.cell[i].invA1[0][1] * sdydU[i].energy;
      grad_E[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].energy +
                    grid.cell[i].invA1[1][1] * sdydU[i].energy;
      
      // Apply angular momentum constraint
      if(param.ang_mom)
      {
         double c[4] = {-grid.cell[i].beta,
                        -grid.cell[i].gamma,
                         grid.cell[i].alpha,
                         grid.cell[i].beta};
         double invAc[4];
         invAc[0] = grid.cell[i].invA1[0][0] * c[0] + grid.cell[i].invA1[0][1] * c[1];
         invAc[1] = grid.cell[i].invA1[1][0] * c[0] + grid.cell[i].invA1[1][1] * c[1];
         invAc[2] = grid.cell[i].invA1[0][0] * c[2] + grid.cell[i].invA1[0][1] * c[3];
         invAc[3] = grid.cell[i].invA1[1][0] * c[2] + grid.cell[i].invA1[1][1] * c[3];
         
         double num = c[0]*grad_rhoU[i].x + c[1]*grad_rhoU[i].y +
                      c[2]*grad_rhoV[i].x + c[3]*grad_rhoV[i].y - ang_mom[i];
         double den = c[0]*invAc[0] + c[1]*invAc[1] + c[2]*invAc[2] + c[3]*invAc[3];
         double lam = num/den;
         
         grad_rhoU[i].x -= lam * invAc[0];
         grad_rhoU[i].y -= lam * invAc[1];
         grad_rhoV[i].x -= lam * invAc[2];
         grad_rhoV[i].y -= lam * invAc[3];
      }
   }
   
   // minmax limiter
   if(param.reconstruct_scheme == Parameter::bj)
      limit_gradients_bj ();
}

//------------------------------------------------------------------------------
// Compute inviscid residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_inviscid_residual ()
{
   // Loop over interior faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      int cl = grid.face[i].lcell;
      int cr = grid.face[i].rcell;
      
      // reconstruct left/right state
      vector<ConVar> state(2);
      reconstruct ( i, grid.face[i].centroid, state );
      
      // If boundary face, apply bc
      if(cr < 0)
      {
         int face_type = grid.face[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         bc.apply (grid.face[i].centroid, grid.face[i], state);
      }
      
      FluxData data;
      //data.ssw = ssw[cl]+ssw[cr];
      //data.ducros = max(ducros[cl], ducros[cr]);

      Flux flux;
      param.material.num_flux (state[0], state[1],
                               grid.face[i].normal, data, flux );

      residual[cl] += flux;
      
      if(cr >= 0)
         residual[cr] -= flux;
      
      if(param.ang_mom)
      {
         Vector drl = grid.face[i].centroid - grid.cell[cl].centroid;
         res_ang_mom[cl] += drl ^ flux.momentum_flux;
         
         if(cr >= 0)
         {
            Vector drr = grid.face[i].centroid - grid.cell[cr].centroid;
            res_ang_mom[cr] -= drr ^ flux.momentum_flux;
         }
      }
   }
}

//------------------------------------------------------------------------------
// Compute inviscid residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_inviscid_residual_2 ()
{
   // Loop over interior faces and accumulate flux
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      int cl = grid.face[i].lcell;
      int cr = grid.face[i].rcell;
      
      unsigned int v0 = grid.face[i].vertex[0];
      unsigned int v1 = grid.face[i].vertex[1];
      
      Vector p0 = grid.vertex[v0].coord * (1.0 - gauss_x[0])/2 + grid.vertex[v1].coord * (1.0 + gauss_x[0])/2;
      Vector p1 = grid.vertex[v0].coord * (1.0 - gauss_x[1])/2 + grid.vertex[v1].coord * (1.0 + gauss_x[1])/2;
      
      // reconstruct left/right state
      vector<ConVar> state0(2), state1(2);
      reconstruct ( i, p0, state0 );
      reconstruct ( i, p1, state1 );
      
      // If boundary face, apply bc
      if(cr < 0)
      {
         int face_type = grid.face[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         bc.apply (p0, grid.face[i], state0);
         bc.apply (p1, grid.face[i], state1);
      }
      
      FluxData data;
      //data.ssw = ssw[cl]+ssw[cr];
      //data.ducros = max(ducros[cl], ducros[cr]);
      
      Flux flux0, flux1;
      param.material.num_flux (state0[0], state0[1],
                               grid.face[i].normal, data, flux0 );
      param.material.num_flux (state1[0], state1[1],
                               grid.face[i].normal, data, flux1 );
      
      flux0 *= gauss_wt[0];
      flux1 *= gauss_wt[1];
      
      residual[cl] += flux0;
      residual[cl] += flux1;
      
      if(cr >= 0)
      {
         residual[cr] -= flux0;
         residual[cr] -= flux1;
      }
      
      if(param.ang_mom)
      {
         Vector drl0 = p0 - grid.cell[cl].centroid;
         res_ang_mom[cl] += drl0 ^ flux0.momentum_flux;
         
         Vector drl1 = p1 - grid.cell[cl].centroid;
         res_ang_mom[cl] += drl1 ^ flux1.momentum_flux;

         if(cr >= 0)
         {
            Vector drr0 = p0 - grid.cell[cr].centroid;
            res_ang_mom[cr] -= drr0 ^ flux0.momentum_flux;
            
            Vector drr1 = p1 - grid.cell[cr].centroid;
            res_ang_mom[cr] -= drr1 ^ flux1.momentum_flux;
         }
      }
   }
}

//------------------------------------------------------------------------------
// Compute residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual ()
{

   // Set residual vector to zero
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      residual[i].zero ();
      res_ang_mom[i] = 0.0;
   }

   // Compute cell and vertex gradients
   compute_gradients ();

   // Inviscid fluxes
   compute_inviscid_residual ();
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FiniteVolume::compute_dt ()
{
   if (param.time_step > 0.0)
   {
      dt_global = param.time_step;
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         dt[i] = dt_global;
      return;
   }

   for(unsigned int i=0; i<grid.n_vertex; ++i)
      dt[i] = 0.0;

   // Interior faces
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      double area = grid.face[i].measure;

      int cl = grid.face[i].lcell;
      double vel_normal_left = primitive[cl].velocity * grid.face[i].normal;
      double c_left = param.material.sound_speed (primitive[cl]);
      dt[cl] += fabs(vel_normal_left) + c_left * area;

      int cr = grid.face[i].rcell;
      if(cr >= 0)
      {
         double vel_normal_right = primitive[cr].velocity * grid.face[i].normal;
         double c_right = param.material.sound_speed (primitive[cr]);
         dt[cr] += fabs(vel_normal_right) + c_right * area;
      }
   }

   // Compute global time step
   dt_global = 1.0e20;

   if( param.time_scheme != "lusgs") 
   {
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         dt[i] = param.cfl * grid.cell[i].area / dt[i];
         dt_global = min( dt_global, dt[i] );
      }
   }
   else
      dt_global = 1.0;

   // For unsteady simulation, use global time step
   if(param.time_mode == "unsteady")
   {
      // Adjust time step so that final time is exactly reached
      if(elapsed_time + dt_global > param.final_time)
         dt_global = param.final_time - elapsed_time;
      for(unsigned int i=0; i<grid.n_cell; ++i)
         dt[i] = dt_global;
   }
}

//------------------------------------------------------------------------------
// Store old conserved variables for multi-stage RK
//------------------------------------------------------------------------------
void FiniteVolume::store_conserved_old ()
{
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      conserved_old[i] = conserved[i];
      ang_mom_old[i] = ang_mom[i];
   }
}

//------------------------------------------------------------------------------
// Update solution to new time level
//------------------------------------------------------------------------------
void FiniteVolume::update_solution (const unsigned int r)
{
   double factor;

   if(param.time_scheme == "ssprk3")
   {
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         factor      = dt[i] / grid.cell[i].area;
         conserved[i]= conserved_old[i] * a_rk[r] +
                       (conserved[i] - residual[i] * factor) * b_rk[r];
         primitive[i]= param.material.con2prim (conserved[i]);
         
         ang_mom[i]= ang_mom_old[i] * a_rk[r] +
                     (ang_mom[i] - res_ang_mom[i] * factor) * b_rk[r];
      }
   }
   else if(param.time_scheme == "rk1" ||
           param.time_scheme == "rk3" ||
           param.time_scheme == "rk4")
   {
      double f = 1.0/(param.n_rks - r);
      for(unsigned int i=0; i<grid.n_cell; ++i)
      {
         factor      = f * dt[i] / grid.cell[i].area;
         conserved[i]= conserved_old[i]  - residual[i] * factor;
         primitive[i]= param.material.con2prim (conserved[i]);
         
         ang_mom[i]= ang_mom_old[i]  - res_ang_mom[i] * factor;
      }
   }
   else if(param.time_scheme == "lusgs")
   { 
      std::cout << "LUSGS not implemented\n";
      exit(0);
   }

}

//------------------------------------------------------------------------------
// Compute L2 norm of mass, momentum and energy residuals
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual_norm (const unsigned int iter)
{
   residual_norm.mass_flux     = 0.0;
   residual_norm.momentum_flux = 0.0;
   residual_norm.energy_flux   = 0.0;

   // Sum of squares for each component
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      double area = grid.cell[i].area;
      residual_norm.mass_flux       += pow(residual[i].mass_flux       / area, 2);
      residual_norm.momentum_flux.x += pow(residual[i].momentum_flux.x / area, 2);
      residual_norm.momentum_flux.y += pow(residual[i].momentum_flux.y / area, 2);
      residual_norm.energy_flux     += pow(residual[i].energy_flux     / area, 2);
   }

   // Take square root and normalize by n_vertex
   residual_norm.mass_flux       = sqrt (residual_norm.mass_flux)       / grid.n_cell;
   residual_norm.momentum_flux.x = sqrt (residual_norm.momentum_flux.x) / grid.n_cell;
   residual_norm.momentum_flux.y = sqrt (residual_norm.momentum_flux.y) / grid.n_cell;
   residual_norm.energy_flux     = sqrt (residual_norm.energy_flux)     / grid.n_cell;

   // Total residual of all components
   residual_norm_total = pow(residual_norm.mass_flux, 2) +
                         residual_norm.momentum_flux.square () +
                         pow(residual_norm.energy_flux, 2);
   residual_norm_total = sqrt (residual_norm_total);

   // Copy residual in first iteration for normalization
   if(iter == last_iter)
   {
      residual_norm_total0 = residual_norm_total;
      cout << "  Initial residual = " << residual_norm_total0 << endl;
      if(residual_norm_total0 == 0.0)
      {
         cout << "  WARNING: Initial residual is zero !!!\n";
         cout << "  WARNING: Setting it to 1.0\n";
         residual_norm_total0 = 1.0;
      }
   }

   residual_norm_total /= residual_norm_total0;
}

//------------------------------------------------------------------------------
// Log messages to screen and file
//------------------------------------------------------------------------------
void FiniteVolume::log_messages (const unsigned int iter)
{

   if(param.time_mode == "steady")
   {
      // File output
      res_file  << setw(8) << iter << "  " 
                << scientific
                << setprecision (4) 
                << dt_global << "  " 
                << residual_norm_total << "  "
                << residual_norm.mass_flux << "  "
                << residual_norm.momentum_flux.x << "  "
                << residual_norm.momentum_flux.y << "  "
                << residual_norm.energy_flux
                << endl;

      // Screen output
      cout << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << residual_norm_total << "  "
           << residual_norm.mass_flux << "  "
           << residual_norm.momentum_flux.x << "  "
           << residual_norm.momentum_flux.y << "  "
           << residual_norm.energy_flux
           << endl;
   }
   else
   {
      // File output
      res_file  << setw(8) << iter << "  " 
                << scientific
                << setprecision (4) 
                << dt_global << "  " 
                << elapsed_time << "  "
                << residual_norm_total
                << endl;

      // Screen output
      cout << setw(8) << iter << "  " 
           << scientific
           << setprecision (4) 
           << dt_global << "  " 
           << elapsed_time << "  "
           << residual_norm_total
           << endl;
   }

   if(bounds)
      compute_bounds (iter);
}

//------------------------------------------------------------------------------
// Save solution to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output (const unsigned int iter, bool write_variables)
{
   static int counter = 0;

   Writer writer (grid, param.material, param.write_format);
   writer.attach_data (primitive);
   writer.attach_gradient (dU, dV, dW);
   if(param.write_variables.size() > 0 && write_variables == true)
      writer.attach_variables (param.write_variables);
   if(param.ducros)
      writer.attach_data (ducros, "ducros");

   writer.output (counter, elapsed_time);
   if(param.time_mode == "unsteady") ++counter;
}

//------------------------------------------------------------------------------
// Save solution to file for restart
//------------------------------------------------------------------------------
void FiniteVolume::output_restart (int iter)
{
   Writer writer (grid);
   writer.attach_data (primitive);
   writer.output_restart (iter);
}

//------------------------------------------------------------------------------
// Find minimum and maximum values in the solution
//------------------------------------------------------------------------------
void FiniteVolume::compute_bounds (const unsigned int iter)
{
   PrimVar prim_min;
   PrimVar prim_max;

   prim_min.density    =  1.0e20;
   prim_min.velocity.x =  1.0e20;
   prim_min.velocity.y =  1.0e20;
   prim_min.pressure   =  1.0e20;

   prim_max.density    = -1.0e20;
   prim_max.velocity.x = -1.0e20;
   prim_max.velocity.y = -1.0e20;
   prim_max.pressure   = -1.0e20;

   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      prim_min.min(primitive[i]);
      prim_max.max(primitive[i]);
   }

   cout << "\t\t Density :"
        << setw(15) << prim_min.density 
        << setw(15) << prim_max.density << endl;
   cout << "\t\t xVelocity   :"
        << setw(15) << prim_min.velocity.x 
        << setw(15) << prim_max.velocity.x << endl;
   cout << "\t\t yVelocity   :"
        << setw(15) << prim_min.velocity.y 
        << setw(15) << prim_max.velocity.y << endl;
   cout << "\t\t Pressure    :"
        << setw(15) << prim_min.pressure 
        << setw(15) << prim_max.pressure << endl;

   if (prim_min.density < 0.0 ||
       prim_min.pressure    < 0.0)
   {
         output (iter, false);
         exit (0);
   }
}

//------------------------------------------------------------------------------
// Compute some quantities like global KE
//------------------------------------------------------------------------------
void FiniteVolume::compute_global (unsigned int iter)
{
   if (!param.has_global) return;

   global_file << iter << "  " << elapsed_time;
   if(param.global_KE)
   {
      double global_KE = 0;
      for(unsigned int i=0; i<grid.n_cell; ++i)
         global_KE += 0.5 * 
                      primitive[i].density *
                      primitive[i].velocity.square() *
                      grid.cell[i].area;
      double global_Enstrophy = 0;
      for(unsigned int i=0; i<grid.n_cell; ++i)
         global_Enstrophy += 0.5 * 
                             pow( dV_cell[i].x - dU_cell[i].y, 2.0) * 
                             grid.cell[i].area;
      global_file << "  " << global_KE << "  " << global_Enstrophy;
   }

   global_file << endl;

}

//------------------------------------------------------------------------------
// Perform time marching iterations
//------------------------------------------------------------------------------
void FiniteVolume::solve ()
{
   unsigned int iter = last_iter;
   elapsed_time = 0.0;
   residual_norm_total = 1.0e20;
   unsigned int last_output_iter = 0;
   
   //output(0); return;
   
   if(param.time_mode == "unsteady")
   {
      compute_gradients ();
      output (0);
   }

   bool found_stop_file = check_for_stop_file ();
   
   while (residual_norm_total > param.min_residue &&
          iter < param.max_iter+last_iter && 
          elapsed_time < param.final_time && !found_stop_file)
   {
      store_conserved_old ();
      compute_dt ();
      //compute_ssw ();
      //compute_ducros ();
      for(unsigned int r=0; r<param.n_rks; ++r)
      {
         compute_residual ();

         //smooth_residual ();
         if(r == param.n_rks-1)
            compute_residual_norm (iter);
         update_solution (r);
      }

      ++iter;
      elapsed_time += dt_global;
      log_messages (iter);

      //compute_forces (iter);
      //compute_global (iter);
      if(iter % param.write_frequency == 0) 
      {
         output (iter);
         last_output_iter = iter;
         found_stop_file = check_for_stop_file ();
      }
   }

   // Save final solution
   if(iter != last_output_iter)
      output (iter);

   if(param.write_restart) output_restart (iter);
}

//------------------------------------------------------------------------------
// This is where the real work starts
//------------------------------------------------------------------------------
void FiniteVolume::run ()
{
   // Read grid from file
   grid.read (param);

   // Set initial condition
   initialize ();
   
   // If -p flag given on command line, then we stop
   if(preprocess)
      return;

   // Measure time taken for solution
   // Store current time
   time_t start_time = time(NULL);

   // Solve the problem
   solve ();

   time_t end_time = time(NULL);
   double time_hours  = difftime(end_time, start_time)/3600.0;
   double time_minutes= difftime(end_time, start_time)/60.0;
   if(time_hours < 1.0)
      cout << "Time taken for computation = " << time_minutes << " minutes\n";
   else
      cout << "Time taken for computation = " << time_hours << " hours\n";
}
