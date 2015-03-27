#ifndef __FV_H__
#define __FV_H__

#include <fstream>
#include <vector>
#include "parameter.h"
#include "material.h"
#include "grid.h"

extern bool preprocess;

// Main class for finite volume scheme
class FiniteVolume
{
   public:
      FiniteVolume (char* file) 
      { 
         param.file = file;
         param.read ();

         // If we are in preprocess mode, then dont open
         // these files. This avoids overwriting existing 
         // files, e.g., when converting to tec or vtk format.
         if(!preprocess)
         {
            res_file.open ("residue.dat");

            if(param.has_global == true)
               global_file.open ("global.dat");
         }
      };
      ~FiniteVolume () 
      {
         res_file.close ();
         if(global_file.is_open())
            global_file.close ();
      };
      void run ();

   private:
      std::ofstream res_file;
      std::ofstream force_file;
      std::ofstream global_file;
      Parameter param;
      Grid      grid;

      std::vector<PrimVar> primitive;
      std::vector<ConVar>  conserved;
      std::vector<ConVar>  conserved_old;
      std::vector<double>  ang_mom, ang_mom_old, res_ang_mom;
      std::vector<Flux>    residual;
      std::vector<Flux>    residual1;
      std::vector<Flux>    residual2;
      std::vector<Vector>  dT, dU, dV, dW, dP;
      std::vector<Vector>  dT_cell, dU_cell, dV_cell, dW_cell, dP_cell;
      std::vector<ConVar>  sdxdU, sdydU;
      std::vector<Vector>  grad_rho, grad_rhoU, grad_rhoV, grad_E;
      std::vector<double>  ducros;
      std::vector<ConVar> phi;
      std::vector<double>  ssw;
      Flux                 residual_norm;
      double               residual_norm_total;
      double               residual_norm_total0;
      std::vector<double>  dt;
      double               dt_global;
      double               elapsed_time;
      int                  last_iter;

      void reconstruct (const unsigned int&      f,
                        const Vector& p,
                        std::vector<ConVar>&    state) const;
      void reconstruct_first (const unsigned int&      f,
                              std::vector<ConVar>&    state) const;
      void reconstruct_second (const unsigned int&      f,
                               const Vector& p,
                               std::vector<ConVar>&    state) const;
      void reconstruct_limited (const unsigned int&      f,
                                std::vector<ConVar>&    state) const;
      PrimVar limited_slope (const ConVar& ul, const ConVar& ur) const;
      void reconstruct_minmod (const unsigned int&      f,
                               std::vector<ConVar>&    state) const;
      PrimVar minmod_slope (const ConVar& ul, const ConVar& ur) const;
      void reconstruct_minmax(const unsigned int&      f,
                              const Vector& p,
                              std::vector<ConVar>&    state) const;
      void prec_thornber(std::vector<PrimVar>& state) const;
      void compute_ssw();
      void compute_ducros();

      void initialize ();
      void interpolate_vertex ();
      void compute_gradients ();
      void limit_gradients_bj ();
      void limit_gradients_mm ();
      void store_conserved_old ();
      void compute_inviscid_residual ();
      void compute_inviscid_residual_2 ();
      void compute_viscous_residual ();
      void compute_axisymmetric_residual ();
      void compute_residual ();
      void smooth_residual ();
      void compute_dt ();
      void compute_residual_norm (const unsigned int iter);
      void log_messages (const unsigned int iter);
      void update_solution (const unsigned int r);
      void solve ();
      void compute_bounds (const unsigned int iter);
      void output (const unsigned int iter, bool write_variables = true);
      void output_restart (int iter);
      void lusgs ();
      void create_force_face_list ();
      void compute_forces (unsigned int iter);
      void compute_global (unsigned int iter);
      
};

#endif
