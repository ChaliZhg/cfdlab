#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <map>
#include "writer.h"

extern Dimension dim;

using namespace std;

//------------------------------------------------------------------------------
// Add primitive variables defined at vertices
//------------------------------------------------------------------------------
void Writer::attach_data (vector<PrimVar>& data)
{
   assert (!has_primitive);
   vertex_primitive = &data;
   has_primitive = true;
}

//------------------------------------------------------------------------------
// Add data defined at vertices
//------------------------------------------------------------------------------
void Writer::attach_data (vector<double>& data, std::string name)
{
   vertex_data.push_back (&data);
   vertex_data_name.push_back (name);
}

//------------------------------------------------------------------------------
// Specify which variables to write
//------------------------------------------------------------------------------
void Writer::attach_variables (const vector<string>& variables)
{
   if(variables.size() > 0)
      assert (has_primitive);

   for(unsigned int i=0; i<variables.size(); ++i)
      if(variables[i]=="mach")
         write_mach = true;
      else if(variables[i]=="vorticity")
         write_vorticity = true;
      else
      {
         cout << "Writer: unknown variable " << variables[i] << endl;
         exit (0);
      }
}

//------------------------------------------------------------------------------
// Add gradient values; currently only velocity gradients added
//------------------------------------------------------------------------------
void Writer::attach_gradient (vector<Vector>& dU_,
                              vector<Vector>& dV_,
                              vector<Vector>& dW_)
{
   assert (!has_gradient);
   dU = &dU_;
   dV = &dV_;
   dW = &dW_;
   has_gradient = true;
}

//------------------------------------------------------------------------------
// Call output function for saving solution to file
//------------------------------------------------------------------------------
void Writer::output (int counter, double elapsed_time)
{
   string filename;
   if     (counter <= 9)    filename = "sol000";
   else if(counter <= 99)   filename = "sol00";
   else if(counter <= 999)  filename = "sol0";
   else if(counter <= 9999) filename = "sol";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
      exit(0);
   }
   stringstream ss;
   ss << counter;
   filename += ss.str();

   if(format == "vtk")
   {
      filename += ".vtk";
      output_vtk (filename);
   }
   else if(format == "tec")
   {
      filename += ".plt";
      output_tec (elapsed_time, filename);
   }
   cout << "Saving solution into file " << filename << endl;
}

//------------------------------------------------------------------------------
// Write data to vtk file
//------------------------------------------------------------------------------
void Writer::output_vtk (string filename)
{
   ofstream vtk;
   vtk.open (filename.c_str());

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "flo3d" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << grid->n_vertex << "  float" << endl;

   for(unsigned int i=0; i<grid->n_vertex; ++i)
      vtk << grid->vertex[i].coord.x << " " 
          << grid->vertex[i].coord.y << " " 
          << 0.0 << endl;

   vtk << "CELLS  " << grid->n_cell << " " << 4 * grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
      vtk << 3 << "  " 
          << grid->cell[i].vertex[0] << " "
          << grid->cell[i].vertex[1] << " "
          << grid->cell[i].vertex[2] << endl;

   vtk << "CELL_TYPES  " << grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
      vtk << 5 << endl;

   // Write vertex data
   if(vertex_data.size() > 0 || has_primitive) 
      vtk << "CELL_DATA  " << grid->n_cell << endl;

   // If vertex primitive data is available, write to file
   if (has_primitive)
   {
      vtk << "SCALARS pressure float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].pressure << endl;

      vtk << "SCALARS density float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].density << endl;

      vtk << "VECTORS velocity float" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_primitive)[i].velocity.x << "  "
             << (*vertex_primitive)[i].velocity.y << "  "
             << 0.0
             << endl;
   }


   // Write mach number
   if(write_mach)
   {
      vtk << "SCALARS mach float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
      {
         double mach = material->Mach ( (*vertex_primitive)[i] );
         vtk << mach << endl;
      }
   }

   // write vorticity
   if(write_vorticity)
   {
      // Check if gradient information is available
      assert(has_gradient);

      vtk << "SCALARS vorticity float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
      {
         double vorticity = (*dV)[i].x - (*dU)[i].y;
         vtk << vorticity << endl;
      }
   }

   // Write vertex data to file
   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      vtk << "SCALARS  " << vertex_data_name[d] << "  float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_cell; ++i)
         vtk << (*vertex_data[d])[i] << endl;
   }

   vtk.close ();
}

//------------------------------------------------------------------------------
// Write data to tecplot file
//------------------------------------------------------------------------------
void Writer::output_tec (double time, string filename)
{
   // Save grid only the first time this function is called
   static bool write_grid = true;
   if(write_grid)
   {
      cout << "Saving grid into grid.plt" << endl;
      ofstream tec;
      tec.open ("grid.plt");
   
      tec << "FILETYPE = GRID" << endl;
      tec << "VARIABLES = \"X\" \"Y\"" << endl;
   
      tec << "ZONE  DATAPACKING=BLOCK, NODES=" << grid->n_vertex 
         << ", ELEMENTS=" << grid->n_cell << ", ZONETYPE=FETRIANGLE" << endl;
   
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << grid->vertex[i].coord.x << endl;
   
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << grid->vertex[i].coord.y << endl;      

      // Triangles
      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << 1+grid->cell[i].vertex[0] << " "
             << 1+grid->cell[i].vertex[1] << " "
             << 1+grid->cell[i].vertex[2] << endl;
   
      write_grid = false;
   }

   // Save solution into different file
   ofstream tec;
   tec.open (filename.c_str());
   
   tec << "FILETYPE = SOLUTION" << endl;
   tec << "VARIABLES = \"P\" \"RHO\" \"U\" \"V\"";
   if(write_mach)
      tec << " \"Mach\"";
   if(write_vorticity)
      tec << " \"Vorticity\"";
   if(dim == axi)
      tec << " \"Vtheta\"";
   for(unsigned int d=0; d<vertex_data.size(); ++d)
      tec << " \"" << vertex_data_name[d] << "\"" << endl;
   tec << endl;
   
   tec << "ZONE STRANDID=1, SOLUTIONTIME=" << time 
       << ", DATAPACKING=BLOCK, NODES=" << grid->n_vertex 
       << ", ELEMENTS=" << grid->n_cell << ", ZONETYPE=FETRIANGLE" << endl;
   
   // If vertex primitive data is available, write to file
   if (has_primitive)
   {
      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << (*vertex_primitive)[i].pressure << endl;

      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << (*vertex_primitive)[i].density << endl;
      
      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << (*vertex_primitive)[i].velocity.x << endl;
      
      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << (*vertex_primitive)[i].velocity.y << endl;

   }
   
   // Write mach number
   if(write_mach)
   {
      for(unsigned int i=0; i<grid->n_cell; ++i)
      {
         double mach = material->Mach ( (*vertex_primitive)[i] );
         tec << mach << endl;
      }
   }
   
   // write vorticity
   if(write_vorticity)
   {
      // Check if gradient information is available
      assert(has_gradient);
      
      for(unsigned int i=0; i<grid->n_cell; ++i)
      {
         double vorticity = (*dV)[i].x - (*dU)[i].y;
         tec << vorticity << endl;
      }
   }

   // Write vertex data to file
   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      for(unsigned int i=0; i<grid->n_cell; ++i)
         tec << (*vertex_data[d])[i] << endl;
   }
   
   tec.close ();
}

//------------------------------------------------------------------------------
// Write solution for restarting
//------------------------------------------------------------------------------
void Writer::output_restart (int iter)
{
   assert (has_primitive);

   cout << "Saving restart file restart.dat\n";

   ofstream fo;
   fo.open ("restart.dat");
   assert (fo.is_open());

   for(unsigned int i=0; i<grid->n_vertex; ++i)
      fo << scientific
         << (*vertex_primitive)[i].density     << "  "
         << (*vertex_primitive)[i].velocity.x  << "  "
         << (*vertex_primitive)[i].velocity.y  << "  "
         << (*vertex_primitive)[i].pressure    << endl;

   fo << iter << endl;
   fo.close ();
}
