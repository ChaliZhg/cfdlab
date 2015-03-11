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
   residual.resize (grid.n_cell);
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
            >> primitive[i].velocity.z
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
// Compute derivatives of velocity and density at grid vertices
//------------------------------------------------------------------------------
void FiniteVolume::compute_gradients ()
{

   std::vector<ConVar> sdxdU(grid.n_cell), sdydU(grid.n_cell);
   double dx, dy;
   
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      sdxdU[i].density = sdydU[i].density = 0.0;
      sdxdU[i].energy  = sdydU[i].energy  = 0.0;
      sdxdU[i].momentum = 0.0;
      sdydU[i].momentum = 0.0;
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
      grad_rho[i].z = 0.0;
      
      grad_rhoU[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].momentum.x +
                       grid.cell[i].invA1[0][1] * sdydU[i].momentum.x;
      grad_rhoU[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].momentum.x +
                       grid.cell[i].invA1[1][1] * sdydU[i].momentum.x;
      grad_rhoU[i].z = 0.0;

      grad_rhoV[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].momentum.y +
                       grid.cell[i].invA1[0][1] * sdydU[i].momentum.y;
      grad_rhoV[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].momentum.y +
                       grid.cell[i].invA1[1][1] * sdydU[i].momentum.y;
      grad_rhoV[i].z = 0.0;

      
      grad_E[i].x = grid.cell[i].invA1[0][0] * sdxdU[i].energy +
                    grid.cell[i].invA1[0][1] * sdydU[i].energy;
      grad_E[i].y = grid.cell[i].invA1[1][0] * sdxdU[i].energy +
                    grid.cell[i].invA1[1][1] * sdydU[i].energy;
      grad_E[i].z = 0.0;
   }

   return;
   
   // TODO
   // minmax limiter
   if(param.reconstruct_scheme == Parameter::bj)
      limit_gradients_bj ();
   else if(param.reconstruct_scheme == Parameter::minmax)
      limit_gradients_mm ();

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
      reconstruct ( i, state );
      
      // If boundary face, apply bc
      if(cr < 0)
      {
         int face_type = grid.face[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         bc.apply (grid.face[i].centroid, grid.face[i], state);
//         std::cout << state[0].momentum.x << " " << state[0].momentum.y << "\n";
//         std::cout << state[1].momentum.x << " " << state[1].momentum.y << "\n";
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
      
//      std::cout << state[0].density << " " << state[0].energy << "\n";
//      std::cout << state[1].density << " " << state[1].energy << "\n";
//      std::cout << grid.face[i].normal.x << "  " << grid.face[i].normal.y << "\n";
//      std::cout << flux.mass_flux << " " << flux.energy_flux << std::endl;
   }
}




//------------------------------------------------------------------------------
// Compute residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual ()
{

   // Set residual vector to zero
   for(unsigned int i=0; i<grid.n_cell; ++i)
      residual[i].zero ();

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
      conserved_old[i] = conserved[i];
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
         //std::cout << residual[i].mass_flux << " " << residual[i].energy_flux << "\n";
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
      }
   }
   else if(param.time_scheme == "lusgs")
   { 
      std::cout << "LUSGS not implemented\n";
      exit(0);
      // Forward Sweep and backward sweep
      lusgs();

      for (unsigned int i=0; i<grid.n_cell; ++i)
      {
         conserved[i] = conserved_old[i] + residual[i];
         primitive[i] = param.material.con2prim (conserved[i]);
      }
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
      residual_norm.momentum_flux.z += pow(residual[i].momentum_flux.z / area, 2);
      residual_norm.energy_flux     += pow(residual[i].energy_flux     / area, 2);
   }

   // Take square root and normalize by n_vertex
   residual_norm.mass_flux       = sqrt (residual_norm.mass_flux)       / grid.n_cell;
   residual_norm.momentum_flux.x = sqrt (residual_norm.momentum_flux.x) / grid.n_cell;
   residual_norm.momentum_flux.y = sqrt (residual_norm.momentum_flux.y) / grid.n_cell;
   residual_norm.momentum_flux.z = sqrt (residual_norm.momentum_flux.z) / grid.n_cell;
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
                << residual_norm.momentum_flux.z << "  "
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
           << residual_norm.momentum_flux.z << "  "
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

   Writer writer (grid, param.material, param.write_format, param.write_surfaces);
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
   prim_min.velocity.z =  1.0e20;
   prim_min.pressure   =  1.0e20;

   prim_max.density    = -1.0e20;
   prim_max.velocity.x = -1.0e20;
   prim_max.velocity.y = -1.0e20;
   prim_max.velocity.z = -1.0e20;
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
   cout << "\t\t zVelocity   :"
        << setw(15) << prim_min.velocity.z 
        << setw(15) << prim_max.velocity.z << endl;
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
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         global_KE += 0.5 * 
                      primitive[i].density *
                      primitive[i].velocity.square() *
                      grid.dcarea[i];
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

   create_force_face_list ();

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
