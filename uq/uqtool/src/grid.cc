#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include "grid.h"

using namespace std;

// Sample constructor
template <int dim>
Sample<dim>::Sample(const unsigned int n_var,
                    const unsigned int n_cell,
                    const unsigned int n_moment,
                    const unsigned int counter)
   :
   n_var  (n_var),
   n_cell (n_cell),
   idx    (counter)
{
   status = NEW;
   load   = false;
   
   J.resize (n_moment);
   
   // Set directory for sample: S000 to S999
   if(counter <= 9)
      sprintf(directory, "S00%d", counter);
   else if(counter <= 99)
      sprintf(directory, "S0%d", counter);
   else if(counter <= 999)
      sprintf(directory, "S%d", counter);
   else
   {
      cout << "Sample: exceeding 1000 !!!" << endl;
      abort ();
   }
}

// Read primal/adjoint solution from file
template <int dim>
void Sample<dim>::read ()
{
   // Make sure solution is not already loaded into memory
   assert (load == false);
   
   primal  = new double [n_var * n_cell];
   adjoint = new double [n_var * n_cell];
   
   // Read from file
   char filename[128];   
   ifstream ff;
   unsigned int c;
   
   // Read primal solution
   sprintf(filename, "RESULT/%s/primal.dat", directory);
   ff.open (filename);
   assert (ff.is_open());
   c = 0;
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         ff >> primal[c++];
   ff.close ();
   
   // Read adjoint solution
   sprintf(filename, "RESULT/%s/adjoint.dat", directory);
   ff.open (filename);
   assert (ff.is_open());
   c = 0;
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         ff >> adjoint[c++];
   ff.close ();
   
   load = true;
}

// Free memory for sample
template <int dim>
void Sample<dim>::clear ()
{
   // Make sure solution is in memory before clearing it
   assert (load == true);
   
   delete [] primal;
   delete [] adjoint;
   
   load = false;
}

// Element constructor
template <int dim>
Element<dim>::Element (const unsigned int order,
                       const unsigned int n_moment,
                       const unsigned int n_cell,
                       const unsigned int counter)
   :
   order    (order),
   n_moment (n_moment),
   n_cell   (n_cell)
{
   assert (order == 1 || order == 2);
   assert (n_moment > 0);
   
   status      = NEW;
   active      = true;
   refine_flag = false;
   moment.resize  (n_moment);
   adj_cor.resize (n_moment);
   RE.resize      (n_moment);

   if(dim == 1)
   {
      if(order==1) n_dof = 2;
      if(order==2) n_dof = 3;
   }
   else if(dim == 2)
   {
      if(order==1) n_dof = 3;
      if(order==2) n_dof = 6;
   }
   else 
   {
      cout << "Element::init : dim = " << dim
           << " not implemented" << endl;
      abort ();
   }
   
   idof.resize(n_dof);
   dof.resize (n_dof);
   
   // Set directory for element: E000 to E999
   if(counter <= 9)
      sprintf(directory, "E00%d", counter);
   else if(counter <= 99)
      sprintf(directory, "E0%d", counter);
   else if(counter <= 999)
      sprintf(directory, "E%d", counter);
   else
   {
      cout << "Element: exceeding 1000 !!!" << endl;
      abort ();
   }
}

// Save mesh_error into files
template <int dim>
void Element<dim>::save_mesh_error ()
{
   unsigned int c = 0;
   
   for(unsigned int i=0; i<n_moment; ++i)
   {
      // Create directory
      char command[128];
      sprintf(command, "mkdir -p RESULT/%s", directory);
      system(command);
      
      // Write mesh error indicator
      char filename[64];
      sprintf(filename, "RESULT/%s/error%d.dat", directory, i);
      ofstream fo;
      fo.open (filename);
      fo.precision (15);
      fo.setf (ios::scientific);
      
      for(unsigned int j=0; j<n_cell; ++j)
         fo << mesh_error[c++] << endl;
      
      fo.close ();
   }
   
   mesh_error.resize (0);
}

// Read element mesh_error from files
template <int dim>
void Element<dim>::load_mesh_error ()
{
   mesh_error.resize (n_moment * n_cell);
   
   unsigned int c = 0;
   
   for(unsigned int i=0; i<n_moment; ++i)
   {
      char filename[64];
      sprintf(filename, "RESULT/%s/error%d.dat", directory, i);
      ifstream fi;
      fi.open (filename);
      assert (fi.is_open());
      
      for(unsigned int j=0; j<n_cell; ++j)
         fi >> mesh_error[c++];
      
      fi.close ();
   }
}

// Set pointer from element.dof to sample
// This must be done every time sample is modified by push_back,
// usually after grid refinement
template <int dim>
void Grid<dim>::reinit_dof (vector<typename Sample<dim>::Sample>& sample)
{
   for(unsigned int i=0; i<element.size(); ++i)
      for(unsigned int j=0; j<element[i].n_dof; ++j)
      {
         unsigned int idx = element[i].idof[j];
         element[i].dof[j] = &sample[idx];
      }
}

// To avoid linker errors
template class Sample<1>;
template class Element<1>;
template class Grid<1>;