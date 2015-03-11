#ifndef __WRITER_H__
#define __WRITER_H__

#include <vector>
#include <string>
#include "grid.h"
#include "material.h"

class Writer
{
   public:
      Writer (const Grid&     grid)
         : 
         grid (&grid),
         has_primitive (false),
         has_gradient (false),
         write_mach (false),
         write_vorticity (false)
         {};
      Writer (const Grid&      grid,
              const Material&  material,
              std::string      format) 
         : 
         grid (&grid),
         material (&material),
         format (format),
         has_primitive (false),
         has_gradient (false),
         write_mach (false),
         write_vorticity (false)
         {};
      void attach_data (std::vector<PrimVar>& data);
      void attach_data (std::vector<double>& data, std::string name);
      void attach_gradient (std::vector<Vector>& dU,
                            std::vector<Vector>& dV,
                            std::vector<Vector>& dW);
      void attach_variables (const std::vector<std::string>& variables);
      void output (int counter, double elapsed_time);
      void output_vtk (std::string filename);
      void output_tec (double time, std::string filename);
      void output_restart (int iter);

   private:

      const Grid*      grid;
      const Material*  material;
      std::string      format;

      std::vector< std::vector<double>* > vertex_data;
      std::vector<std::string> vertex_data_name;

      std::vector<PrimVar>* vertex_primitive;
      std::vector<Vector>* dU;
      std::vector<Vector>* dV;
      std::vector<Vector>* dW;
      bool has_primitive;
      bool has_gradient;

      bool write_mach;
      bool write_vorticity;

};

#endif
