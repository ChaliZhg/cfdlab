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
   
   // Read from file TBD
   char filename[128];   
   ifstream ff;
   unsigned int c;
   
   // Read primal solution
   sprintf(filename, "%s/primal.dat", directory);
   ff.open (filename);
   assert (ff.is_open());
   c = 0;
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         ff >> primal[c++];
   ff.close ();
   
   // Read adjoint solution
   sprintf(filename, "%s/adjoint.dat", directory);
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
                       const unsigned int n_moment)
   :
   order    (order),
   n_moment (n_moment)
{
   assert (order == 1 || order == 2);
   assert (n_moment > 0);
   
   status      = NEW;
   active      = true;
   refine_flag = false;
   parent      = this;
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
   }
   
   dof.resize (n_dof);
}

// To avoid linker errors
template class Sample<1>;
template class Element<1>;
template class Grid<1>;