#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include "uq.h"
#include "quadrature.h"
#include "interpolate.h"
#include "evaluator.h"

using namespace std;

// Constructor
// Read options from file uq.in
// Copy initial template directory
template <int dim>
UQProblem<dim>::UQProblem ()
{
   read_options ();
   
   assert (n_moment == 1);
   
   moment.resize  (n_moment);
   adj_cor.resize (n_moment);
   RE.resize      (n_moment);
   
   // Copy initial template directory
   sprintf(template_dir, "template_dir_0");
   char command[128];
   sprintf(command, "cp -r template_dir RESULT/%s", template_dir);
   system (command);
}

// Destructor
template <int dim>
UQProblem<dim>::~UQProblem ()
{
}

// Read options from file uq.in
template <int dim>
void UQProblem<dim>::read_options ()
{   
   
   cout << "Reading input from uq.in ...\n";
   
   ifstream fi;
   fi.open ("uq.in");
   assert (fi.is_open());
   
   int dim_in;
   string str;
   
   fi >> str >> dim_in;
   assert (str == "dim");
   assert (dim == dim_in);
   
   string input;
   
   for(unsigned int i=0; i<dim; ++i)
   {
      fi >> pdf_data.x_name[i]
         >> pdf_data.x_min[i]
         >> pdf_data.x_max[i]
         >> input;
      if(input=="uniform") 
         pdf_data.type[i] = uniform;
      else if(input=="normal")  
      {
         pdf_data.type[i] = normal;
         // Mean value is at center of given range
         pdf_data.mean[i] = 0.5 * (pdf_data.x_min[i] + pdf_data.x_max[i]); 
         fi >> input >> pdf_data.variance[i];
         assert (input == "variance");
         assert (pdf_data.variance[i] > 0.0);
      }
      else
      {
         cout << "Unknown PDF : " << input << endl;
         abort ();
      }
   }
   
   fi >> str >> n_moment;
   assert (str == "n_moment");
   assert (n_moment >= 1);
   
   fi >> str >> n_var;
   assert (str == "n_var");
   assert (n_var >= 1);
   
   fi >> str >> n_cell;
   assert (str == "n_cell");
   assert (n_cell >= 1);
   
   fi >> str >> max_sample;
   assert (str == "max_sample");
   assert (max_sample > 0);
   
   // Error control type
   fi >> str >> input;
   assert (str == "error");
   if(input=="stochastic")
      error_control = STOCHASTIC;
   else if(input=="combined")
      error_control = COMBINED;
   else
   {
      cout << "Unknown error control type : " << input << endl;
      abort ();
   }
   
   // Refinement type
   fi >> str >> input;
   assert (str == "refine");
   if(input=="uniform")
      refine_type = UNIFORM;
   else if(input=="adaptive")
      refine_type = ADAPTIVE;
   else
   {
      cout << "Unknown refinement type : " << input << endl;
      abort ();
   }
   
   fi.close ();
}

// Perform new primal/adjoint simulations
template <int dim>
void UQProblem<dim>::run_simulations ()
{
   for(unsigned int i=0; i<sample.size(); ++i)
      if(sample[i].status == NEW)
      {         
         cout << "Running simulation for sample " << sample[i].directory << endl;
         
         // Create directory if it does not exist
         char command[256];
         sprintf(command, "mkdir -p RESULT/%s", sample[i].directory);
         system (command);

         // Print random variables to file
         char filename[64];
         sprintf(filename, "RESULT/%s/random.dat", sample[i].directory);
         ofstream ff;
         ff.open (filename);
         ff.precision (15);
         ff.setf (ios::scientific);
         for(unsigned int j=0; j<dim; ++j)
            ff << "{" << pdf_data.x_name[j] << "}  " <<
                  sample[i].x[j] << endl;
         ff.close ();
         
         // Run external code inside the sample directory
         // TBD Check that simulation was successful
         sprintf(command, "./runsolver.sh 1 RESULT/%s RESULT/%s", 
                 sample[i].directory, template_dir);
         system (command);
         
         // Read objective functions
         ifstream fi;
         sprintf(filename, "RESULT/%s/obj.dat", sample[i].directory);
         fi.open (filename);
         assert (fi.is_open());
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
   JREvaluator<dim> evaluator (pdf_data.x_name, n_moment, n_cell,
                               template_dir);
   
   // Compute moments for newly created elements
   for(unsigned int i=0; i<grid.element.size(); ++i)
      if(grid.element[i].status == NEW)
      {
         cout << "Computing mean for element = " << i << endl;
         
         // Allocate memory for mesh error indicator and set to zero
         if(error_control == COMBINED)
            grid.element[i].mesh_error.resize (n_moment * n_cell, 0.0);
         
         grid.element[i].moment  = 0.0;
         grid.element[i].adj_cor = 0.0;
         grid.element[i].RE      = 0.0;
         
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
            
            grid.element[i].moment  += w * evaluator.J;
            grid.element[i].adj_cor += w * evaluator.VdotR;
            grid.element[i].RE      += w * evaluator.RE;
            
            if(error_control == COMBINED)
               grid.element[i].mesh_error += w * evaluator.RE_array;
            
         }
         
         // Clear all samples belonging to this element from memory
         for(unsigned int d=0; d<grid.element[i].n_dof; ++d)
            grid.element[i].dof[d]->clear();
                     
         // Save mesh_error into file
         if(error_control == COMBINED)
            grid.element[i].save_mesh_error ();

         grid.element[i].status = OLD;
      }

   // Accumulate moment contribution from all active elements   
   moment = adj_cor = RE = 0.0;
   for(unsigned int j=0; j<grid.element.size(); ++j)
      if (grid.element[j].active)
      {
         moment  += grid.element[j].moment;
         adj_cor += grid.element[j].adj_cor;
         for(unsigned int i=0; i<n_moment; ++i)
            RE[i] += fabs(grid.element[j].RE[i]);
      }
   
   // Accumulate physical cell error by summing over all
   // stochastic elements
   if(error_control == COMBINED)
   {
      mesh_error.resize (n_moment * n_cell, 0.0);
      for(unsigned int i=0; i<grid.element.size(); ++i)
         if (grid.element[i].active)
         {
            // Load element mesh_error from file
            grid.element[i].load_mesh_error ();
            mesh_error += grid.element[i].mesh_error;
            grid.element[i].mesh_error.resize(0);
         }
   }

}

// Find element with largest remaining error and flag it for
// stochastic refinement
template <int dim>
void UQProblem<dim>::flag_elements ()
{
   // Uniform refinement: flag all active elements
   if(refine_type == UNIFORM)
   {
      for(unsigned int j=0; j<grid.element.size(); ++j)
         if(grid.element[j].active)
            grid.element[j].refine_flag = true;
   }
   else
   {  // Adaptive refinement
      for(unsigned int i=0; i<n_moment; ++i)
      {
         unsigned int imax = 0;
         double max_error = -1.0;
         
         // Find stochastic element with largest |RE|
         // Only check active elements
         for(unsigned int j=0; j<grid.element.size(); ++j)
            if(grid.element[j].active &&
               fabs(grid.element[j].RE[i]) > max_error)
            {
               imax = j;
               max_error = fabs(grid.element[j].RE[i]);
            }
         
         grid.element[imax].refine_flag = true;
      }
   }
}

// Refine physical grid
template <int dim>
void UQProblem<dim>::refine_physical (const unsigned int iter)
{   
   // Create name for new_template_dir
   char new_template_dir[64];
   sprintf(new_template_dir, "template_dir_%d", iter);
   
   // Copy contents of latest template_dir into new_template_dir
   char command[128];
   sprintf(command, "cp -r RESULT/%s RESULT/%s", template_dir, 
           new_template_dir);
   system (command);
   
   // Save mesh_error into file
   // Separate file is created for each moment
   unsigned int c = 0;
   for(unsigned int i=0; i<n_moment; ++i)
   {
      char filename[64];
      sprintf(filename, "RESULT/%s/error%d.dat", new_template_dir, i);
      ofstream fo;
      fo.open (filename);
      fo.precision (15);
      fo.setf (ios::scientific);
      
      for(unsigned int j=0; j<n_cell; ++j)
         fo << mesh_error[c++] << endl;
      
      fo.close ();
   }
   
   // Call external grid refinement program
   sprintf (command, "./adapt.sh RESULT/%s", new_template_dir);
   system (command);
   
   // Read new number of cells: n_cell
   char filename[64];
   unsigned int n_cell_old, n_cell_new;
   ifstream fi;
   sprintf(filename, "RESULT/%s/cell_map.dat", new_template_dir);
   fi.open (filename);
   assert (fi.is_open());
   fi >> n_cell_old >> n_cell_new;
   fi.close ();
   
   // Just a check
   assert (n_cell_old == n_cell_old);
   
   // If grid has been adapted, set sample status to re-run
   if(n_cell != n_cell_new)
   {
      cout << "Physical mesh has been adapted\n";
      
      n_cell = n_cell_new;
      
      // All samples must be re-run on new grid
      for(unsigned int i=0; i<sample.size(); ++i)
      {
         sample[i].n_cell = n_cell;
         sample[i].status = NEW;
      }
      
      // Moment contribution from all active stochastic elements 
      // must be re-computed
      for(unsigned int i=0; i<grid.element.size(); ++i)
         if(grid.element[i].active)
         {
            grid.element[i].n_cell = n_cell;
            grid.element[i].status = NEW;
         }
   }
   else
   {
      cout << "Physical mesh was not adapted\n";
   }
   
   // Update template dir
   sprintf(template_dir, "%s", new_template_dir);
}

// Print messages to screen
// Write moments etc. to file
template <int dim>
void UQProblem<dim>::log_result (ofstream& fj, ofstream& fe)
{
   // Count number of active elements
   unsigned int n_elem_active = 0;
   for(unsigned int i=0; i<grid.element.size(); ++i)
      if(grid.element[i].active)
         ++n_elem_active;
   
   cout << "No. of stochastic samples  = " << sample.size() << endl;
   cout << "No. of stochastic elements = " << grid.element.size()
        << endl;
   cout << "No. of active elements     = " << n_elem_active << endl;
   cout << "No. of physical cells      = " << n_cell << endl;  
   
   fj.precision(15);
   fj.setf (ios::scientific);
   
   fe.precision(15);
   fe.setf (ios::scientific);
   
   fj << sample.size() << " " << n_elem_active << " ";
   fe << sample.size() << " " 
      << 1.0/pow(double(sample.size()), 1.0/dim) << " ";

   for(unsigned int i=0; i<n_moment; ++i)
   {
      fj << moment[i] << " "
         << moment[i] + adj_cor[i] << " ";
      fe << fabs(adj_cor[i]) << " "
         << fabs(RE[i]) << " ";
   }      
   
   fj << endl; fj.flush ();
   fe << endl; fe.flush ();

}

// This is where it all begins
template <int dim>
void UQProblem<dim>::run ()
{
   ofstream fj, fe;
   fj.open ("RESULT/J.dat");
   fe.open ("RESULT/error.dat");
   
   make_grid ();

   unsigned int iter = 0;
   
   while (sample.size() <= max_sample)
   {      
      cout << "------------------" << endl;
      cout << "Iteration = " << iter << endl;
      cout << "------------------" << endl;
      
      if (iter > 0) 
      {
         flag_elements ();
         refine_grid ();
         grid.reinit_dof (sample);
         if(error_control == COMBINED)
            refine_physical (iter);
      }
      run_simulations ();
      compute_moments ();

      log_result (fj, fe);   
      output (iter);
      ++iter;
   }
   
   fj.close ();
   fe.close ();

}



// To avoid linker errors
template class UQProblem<1>;
