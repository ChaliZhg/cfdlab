#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include "uq.h"
#include "quadrature.h"
#include "interpolate.h"
#include "evaluator.h"

using namespace std;

// Constructor
template <int dim>
UQProblem<dim>::UQProblem ()
{
   read_options ();
   
   moment  = new double [n_moment];
   adj_cor = new double [n_moment];
}

// Destructor
template <int dim>
UQProblem<dim>::~UQProblem ()
{
   delete [] moment;
   delete [] adj_cor;
}

// Read options from file
// For now, we set it here itself
template <int dim>
void UQProblem<dim>::read_options ()
{   
   n_moment = 1;
   
   pdf_data.x_min[0] = 0.0;
   pdf_data.x_max[0] = 1.0;
   pdf_data.type[0]  = uniform;
   
   // Determinitic simulation size
   n_var = 1;
   n_cell= 50;
   
   max_sample = 100;
}

// Perform new primal/adjoint simulations
template <int dim>
void UQProblem<dim>::run_simulations ()
{
   for(unsigned int i=0; i<sample.size(); ++i)
      if(sample[i].status == NEW)
      {         
         // Create directory if it does not exist
         char command[128];
         sprintf(command, "mkdir -f %s", sample[i].directory);
         system (command);

         // Print random variables to file
         char filename[64];
         sprintf(filename, "%s/random.dat", sample[i].directory);
         ofstream ff;
         ff.open (filename);
         for(unsigned int j=0; j<dim; ++j)
            ff << sample[i].x[j] << endl;
         ff.close ();
         
         // Run external code inside the sample directory
         // TBD Check that simulation was successful
         sprintf(command, "./runsolver.sh %s", sample[i].directory);
         system (command);
         
         // Read objective functions
         ifstream fi;
         sprintf(filename, "%s/obj.dat", sample[i].directory);
         fi.open (filename);
         for(unsigned int j=0; j<n_moment; ++j)
            fi >> sample[i].J[j];
         fi.close ();
         
         sample[i].status = OLD;
      }
}

// Compute moments and adjoint correction on each new stochastic element
template <int dim>
void UQProblem<dim>::compute_moments ()
{
   Quadrature<dim>  quadrature_formula;
   Interpolate<dim> interpolate_formula (n_var, n_cell);
   JREvaluator<dim> evaluator (n_moment);
   
   // Compute moments for newly created elements
   for(unsigned int i=0; i<grid.n_element; ++i)
      if(grid.element[i].status == NEW)
      {
         for(unsigned int m=0; m<n_moment; ++m)
         {
            grid.element[i].moment[m]  = 0.0;
            grid.element[i].adj_cor[m] = 0.0;
         }
         
         // Read all samples belonging to this element from file
         for(unsigned int d=0; d<grid.element[i].n_dof; ++d)
            grid.element[i].dof[d]->read();
         
         // Perform quadrature
         quadrature_formula.reinit(grid.element[i]);
         interpolate_formula.reinit(grid.element[i]);
         unsigned int n_q_point = quadrature_formula.n_point;
         for(unsigned int q=0; q<n_q_point; ++q)
         {
            double* x = quadrature_formula.q_point[q];
            interpolate_formula.execute (x);
            // Compute integrand, primal/adjoint residual
            evaluator.execute (x, interpolate_formula);
            double pdf = pdf_data.get_pdf (x);
            double w   = pdf * quadrature_formula.weight[q];
            
            for(unsigned int m=0; m<n_moment; ++m)
            {
               grid.element[i].moment[m]  += w * evaluator.J[m];
               grid.element[i].adj_cor[m] += w * evaluator.VdotR[m];

            }
         }
         
         // Clear all samples belonging to this element from memory
         for(unsigned int d=0; d<grid.element[i].n_dof; ++d)
            grid.element[i].dof[d]->clear();

         grid.element[i].status = OLD;
      }

   // Accumulate moment contribution from all active elements
   for(unsigned int i=0; i<n_moment; ++i)
   {
      moment[i] = adj_cor[i] = 0.0;
      for(unsigned int j=0; j<grid.n_element; ++j)
         if (grid.element[j].active)
         {
            moment[i]  += grid.element[j].moment[i];
            adj_cor[i] += grid.element[j].adj_cor[i];
         }
   }
}

// This is where it all begins
template <int dim>
void UQProblem<dim>::run ()
{
   make_grid ();

   unsigned int iter = 0;
   
   while (n_sample < max_sample)
   {
      if (iter > 0) refine_grid ();
      run_simulations ();
      compute_moments ();
      
      ++iter;

      cout << "Iteration = " << iter << endl;
      cout << "No. of samples = " << n_sample << endl;
      cout << "No. of cells   = " << n_cell << endl;
   }

}

// To avoid linker errors
template class UQProblem<1>;
