#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cassert>
#include "grid.h"
#include "parameter.h"

using namespace std;

//------------------------------------------------------------------------------
// Read grid from file
//------------------------------------------------------------------------------
void Grid::read (const Parameter& param)
{
   if(param.grid_type == gmsh)
      read_gmsh (param.grid_file);
   else if(param.grid_type == bamg)
      read_bamg (param.grid_file);
   else if(param.grid_type == delaundo)
      read_delaundo (param.grid_file);
   else
   {
      cout << "Unknown grid type specified !!!" << endl;
      exit (0);
   }

   /*
   int c = 0;
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      double yy = fabs(vertex[i].coord.y - 0.5);
      if(yy < 1.0e-10)
      {
         vertex[i].coord.y += pow(-1.0, c) * 1e-3;
         ++c;
      }
   }
   printf("Number of perturbed points = %d\n", c);
   */

   // At this stage, we have only boundary faces. We save this number.
   n_boundary_face = n_face;

   check_face_type (param.boundary_condition);
   preproc ();
   info ();
}

//------------------------------------------------------------------------------
// Print some grid information to screen
//------------------------------------------------------------------------------
void Grid::info ()
{

   double min_face_length =  1.0e20;
   double max_face_length = -1.0e20;

   for(unsigned int i=0; i<n_face; ++i)
   {
      min_face_length = min ( min_face_length, face[i].normal.norm() );
      max_face_length = max ( max_face_length, face[i].normal.norm() );
   }
   
   double min_cell_area =  1.0e20;
   double max_cell_area = -1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      min_cell_area = min ( min_cell_area, cell[i].area );
      max_cell_area = max ( max_cell_area, cell[i].area );
   }

   cout << "Grid information:\n";
   cout << "   Number of vertices       = " << n_vertex << endl;
   cout << "   Number of triangles      = " << n_cell << endl;
   cout << "   Number of boundary edges = " << n_boundary_face << endl;
   cout << setw(30) << "min" << setw(15) << "max" << endl;
   cout << "  cell area    =  " << setw(15) << min_cell_area 
                                << setw(15) << max_cell_area << endl;
   cout << "  face length  =  " << setw(15) << min_face_length
                                << setw(15) << max_face_length << endl;
}

//------------------------------------------------------------------------------
// Check that all boundary faces have been assigned a bc type
//------------------------------------------------------------------------------
void Grid::check_face_type (const map<int,BoundaryCondition>& bc)
{
   for(unsigned int i=0; i<face.size(); ++i)
   {
      assert (face[i].type != -1);
      if(bc.find(face[i].type) == bc.end())
      {
         cout << "check_face_type:\n";
         cout << "   No boundary condition specified for\n";
         cout << "   face = " << i << " whose type = " << face[i].type << endl;
         cout << "   There may be more faces with similar problem.\n";
         exit (0);
      }
   }
}
