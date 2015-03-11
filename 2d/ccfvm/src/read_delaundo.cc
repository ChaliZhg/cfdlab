#include <iostream>
#include <fstream>
#include <cassert>
#include <string.h>
#include <cstdlib>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Read a grid in delaundo format
//------------------------------------------------------------------------------
void Grid::read_delaundo (string grid_file)
{
   unsigned int i, idummy;
   double rdummy;
   string input;

   cout << "Reading delaundo grid file " << grid_file << endl;

   ifstream file;
   file.open (grid_file.c_str());
   assert ( file.is_open() );

   getline(file, input);

   // Read triangles
   file >> n_cell
        >> idummy 
        >> idummy;
   assert(n_cell > 0);
   cell.resize(n_cell);
   for(i=0;i<n_cell;i++)
   {
      file >> idummy
           >> cell[i].vertex[0] 
           >> cell[i].vertex[1] 
           >> cell[i].vertex[2] 
           >> idummy
           >> idummy
           >> idummy
           >> idummy;
      --cell[i].vertex[0];
      --cell[i].vertex[1];
      --cell[i].vertex[2];
   }

   // Read vertices
   file >> n_vertex;
   assert(n_vertex > 0);
   vertex.resize(n_vertex);

   // skip a line
   file >> rdummy 
        >> rdummy 
        >> rdummy 
        >> rdummy 
        >> rdummy 
        >> rdummy;

   for(i=0; i<n_vertex; i++)
   {
      file >> vertex[i].coord.x 
           >> vertex[i].coord.y 
           >> rdummy
           >> rdummy
           >> rdummy
           >> rdummy
           >> idummy;
   }

   unsigned int n_seg, n_corner, c=0;
   file >> n_seg >> n_corner;
   for(unsigned int s=0; s<n_corner; ++s)
      file >> idummy;
   for(unsigned int s=0; s<n_corner; ++s)
      file >> idummy;
   for(unsigned int s=0; s<n_seg; ++s)
   {
      unsigned int ne, etype;
      file >> ne >> etype;
      assert (ne > 0);
      face.resize( face.size() + ne );
      for(i=0; i<ne; ++i)
      {
         file >> face[c].vertex[0]
              >> face[c].vertex[1]
              >> idummy
              >> idummy;
         --face[c].vertex[0];
         --face[c].vertex[1];
         face[c].type = etype;
         ++c;
      }
   }

   // total number of face
   n_face = face.size();

   file.close();
}
