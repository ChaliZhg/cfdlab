#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "reader.h"
#include "vec.h"
#include "material.h"
#include "ic.h"
#include "bc.h"
#include "constants.h"

class Parameter
{
   public:
      char* file;

      std::string time_mode;
      std::string time_scheme;
      unsigned int n_rks;
      unsigned int max_iter;
      double time_step;
      double cfl;
      double final_time;
      double min_residue;
      bool smooth_res;
      bool ducros;
      bool ang_mom;

      enum ReconstructScheme 
      { 
         first, second, bj
      };

      ReconstructScheme reconstruct_scheme;

      Material material;

      std::string grid_file;
      GridType    grid_type;

      InitialCondition initial_condition;

      std::map<int,BoundaryCondition> boundary_condition;

      std::string  write_format;
      unsigned int write_frequency;
      std::vector<std::string> write_variables;
      bool write_restart;
      bool has_global;
      bool global_KE;

      void read ();

   private:
      void read_constants (Reader&);
      void read_grid (Reader&);
      void read_numeric (Reader&);
      void read_material (Reader&);
      void read_initial_condition (Reader&);
      void read_boundary (Reader&);
      void read_integrals (Reader&);
      void read_output (Reader&);
      void check ();
};

#endif
