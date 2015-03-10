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
      else if(variables[i]=="density")
         write_density = true;
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

   // save solution at specified surfaces
   string surffilename;
   if     (counter <= 9)    surffilename = "000";
   else if(counter <= 99)   surffilename = "00";
   else if(counter <= 999)  surffilename = "0";
   else if(counter <= 9999) surffilename = "";
   else
   {
      cout << "Writer::output: counter is too large !!!\n";
      exit(0);
   }
   surffilename += ss.str();
   output_surfaces (surffilename);
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
          << grid->vertex[i].coord.z << endl;

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
      vtk << "POINT_DATA  " << grid->n_vertex << endl;

   // If vertex primitive data is available, write to file
   if (has_primitive)
   {
      vtk << "SCALARS pressure float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk << (*vertex_primitive)[i].pressure << endl;

      vtk << "SCALARS temperature float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk << (*vertex_primitive)[i].temperature << endl;

      if(dim == axi)
      {
         vtk << "SCALARS Vtheta float 1" << endl;
         vtk << "LOOKUP_TABLE default" << endl;
         for(unsigned int i=0; i<grid->n_vertex; ++i)
            vtk << (*vertex_primitive)[i].velocity.z << endl;
      }

      vtk << "VECTORS velocity float" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
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
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double mach = material->Mach ( (*vertex_primitive)[i] );
         vtk << mach << endl;
      }
   }

   // Write density
   if(write_density)
   {
      vtk << "SCALARS density float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double density = material->Density ((*vertex_primitive)[i]);
         vtk << density << endl;
      }
   }

   // write vorticity
   if(write_vorticity)
   {
      // Check if gradient information is available
      assert(has_gradient);

      vtk << "SCALARS vorticity float 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
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
      for(unsigned int i=0; i<grid->n_vertex; ++i)
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
   tec << "VARIABLES = \"P\" \"T\" \"U\" \"V\"";
   if(write_mach)
      tec << " \"Mach\"";
   if(write_density)
      tec << " \"Density\"";
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
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_primitive)[i].pressure << endl;

      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_primitive)[i].temperature << endl;
      
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_primitive)[i].velocity.x << endl;
      
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_primitive)[i].velocity.y << endl;

   }
   
   
   // Write mach number
   if(write_mach)
   {
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double mach = material->Mach ( (*vertex_primitive)[i] );
         tec << mach << endl;
      }
   }
   
   // Write density
   if(write_density)
   {
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double density = material->Density ((*vertex_primitive)[i]);
         tec << density << endl;
      }
   }
   
   // write vorticity
   if(write_vorticity)
   {
      // Check if gradient information is available
      assert(has_gradient);
      
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double vorticity = (*dV)[i].x - (*dU)[i].y;
         tec << vorticity << endl;
      }
   }

   // Theta component of velocity
   if(dim == axi)
   {
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_primitive)[i].velocity.z << endl;
   }

   // Write vertex data to file
   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         tec << (*vertex_data[d])[i] << endl;
   }
   
   tec.close ();
}

//------------------------------------------------------------------------------
// Write solution at surfaces
// For each type, two files are created, one with data at vertices (pressure)
// and another at face centers (skin friction)
//------------------------------------------------------------------------------
void Writer::output_surfaces (string surffilename)
{
   if(surfaces.size() == 0) return;

   const int nsurf = surfaces.size();
   vector<ofstream*> fv (nsurf);
   vector<ofstream*> ff (nsurf);
   map<int,int> type_to_idx;
   for(int i=0; i<nsurf; ++i)
   {
      stringstream ss;
      ss << surfaces[i];
      string vfilename = "v" + surffilename + "_" + ss.str() + ".dat";
      string ffilename = "f" + surffilename + "_" + ss.str() + ".dat";
      fv[i] = new ofstream(vfilename.c_str());
      ff[i] = new ofstream(ffilename.c_str());
      assert (fv[i]->is_open());
      assert (ff[i]->is_open());
      type_to_idx.insert(pair<int,int>(surfaces[i], i));
   }
   for(unsigned int i=0; i<grid->bface.size(); ++i)
   {
      const int type = grid->bface[i].type;
      const unsigned int v0 = grid->bface[i].vertex[0];
      const unsigned int v1 = grid->bface[i].vertex[1];
      map<int,int>::iterator it;
      it = type_to_idx.find(type);
      if(it != type_to_idx.end())
      {
         // viscous force, using vertex gradients
         const double T = ((*vertex_primitive)[v0].temperature +
                           (*vertex_primitive)[v1].temperature)/2.0;
         const double mu = material->viscosity(T);
         const Vector gradU = ((*dU)[v0] + (*dU)[v1])/2.0;
         const Vector gradV = ((*dV)[v0] + (*dV)[v1])/2.0;
         const double div = gradU.x + gradV.y;
         const double sxx = 2.0 * mu * (gradU.x - (1.0/3.0) * div);
         const double syy = 2.0 * mu * (gradV.y - (1.0/3.0) * div);
         const double sxy = mu * (gradU.y + gradV.x);
         const Vector normal = grid->bface[i].normal / grid->bface[i].measure;
         const double fx = (sxx * normal.x + sxy * normal.y);
         const double fy = (sxy * normal.x + syy * normal.y);
         const double cf = fx*normal.y - fy*normal.x;
         const Vector center = (grid->vertex[v0].coord + grid->vertex[v1].coord)/2.0;

         const int j = it->second;
         *fv[j] << grid->vertex[v0].coord.x << "  "
                << grid->vertex[v0].coord.y << "  "
                << (*vertex_primitive)[v0].pressure << "  "
                << (*vertex_primitive)[v0].temperature << "  "
                << (*vertex_primitive)[i].velocity.x  << "  "
                << (*vertex_primitive)[i].velocity.y  << "  "
                << (*vertex_primitive)[i].velocity.z  << "  "
                << endl;
         *ff[j] << center.x << "  " << center.y << "  "
                << cf << "  "
                << endl;
      }
   }
   for(int i=0; i<nsurf; ++i)
   {
      fv[i]->close();
      ff[i]->close();
   }
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
         << (*vertex_primitive)[i].temperature << "  "
         << (*vertex_primitive)[i].velocity.x  << "  "
         << (*vertex_primitive)[i].velocity.y  << "  "
         << (*vertex_primitive)[i].velocity.z  << "  "
         << (*vertex_primitive)[i].pressure    << endl;

   fo << iter << endl;
   fo.close ();
}
