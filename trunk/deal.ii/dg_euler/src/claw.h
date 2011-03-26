#ifndef __CLAW_H__
#define __CLAW_H__

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>

#include <lac/vector.h>
#include <lac/compressed_sparsity_pattern.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_in.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>

#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <numerics/solution_transfer.h>

#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_precondition.h>
#include <lac/trilinos_solver.h>


#include <Sacado.hpp>


#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "parameters.h"
#include "integrator.h"

				 // @sect3{Conservation law class}

				 // Here finally comes the class that
				 // actually does something with all
				 // the Euler equation and parameter
				 // specifics we've defined above. The
				 // public interface is pretty much
				 // the same as always (the
				 // constructor now takes the name of
				 // a file from which to read
				 // parameters, which is passed on the
				 // command line). The private
				 // function interface is also pretty
				 // similar to the usual arrangement,
				 // with the
				 // <code>assemble_system</code>
				 // function split into three parts:
				 // one that contains the main loop
				 // over all cells and that then calls
				 // the other two for integrals over
				 // cells and faces, respectively.
template <int dim>
class ConservationLaw
{
public:
   ConservationLaw (const char *input_filename);
   void run ();
   
private:
   void setup_system ();
   
   void setup_mesh_worker (Integrator<dim>&);
   
   std::pair<unsigned int, double> solve (dealii::Vector<double> &solution);
   
   void compute_refinement_indicators (dealii::Vector<double> &indicator) const;
   void refine_grid (const dealii::Vector<double> &indicator);
   
   void output_results () const;
   
   typedef dealii::MeshWorker::DoFInfo<dim> DoFInfo;
   typedef dealii::MeshWorker::IntegrationInfo<dim> CellInfo;
                     
   void integrate_cell_term (DoFInfo& dinfo, CellInfo& info);
   void integrate_boundary_term (DoFInfo& dinfo, CellInfo& info);
   void integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                             CellInfo& info1, CellInfo& info2);
   
   void assemble_system (Integrator<dim>& integrator);
   
   // The first few member variables
   // are also rather standard. Note
   // that we define a mapping
   // object to be used throughout
   // the program when assembling
   // terms (we will hand it to
   // every FEValues and
   // FEFaceValues object); the
   // mapping we use is just the
   // standard $Q_1$ mapping --
   // nothing fancy, in other words
   // -- but declaring one here and
   // using it throughout the
   // program will make it simpler
   // later on to change it if that
   // should become necessary. This
   // is, in fact, rather pertinent:
   // it is known that for
   // transsonic simulations with
   // the Euler equations,
   // computations do not converge
   // even as $h\rightarrow 0$ if
   // the boundary approximation is
   // not of sufficiently high
   // order.
   dealii::Triangulation<dim>   triangulation;
   const dealii::MappingQ1<dim> mapping;
   
   const dealii::FESystem<dim>  fe;
   dealii::DoFHandler<dim>      dof_handler;
   
   // Next come a number of data
   // vectors that correspond to the
   // solution of the previous time
   // step
   // (<code>old_solution</code>),
   // the best guess of the current
   // solution
   // (<code>current_solution</code>;
   // we say <i>guess</i> because
   // the Newton iteration to
   // compute it may not have
   // converged yet, whereas
   // <code>old_solution</code>
   // refers to the fully converged
   // final result of the previous
   // time step), and a predictor
   // for the solution at the next
   // time step, computed by
   // extrapolating the current and
   // previous solution one time
   // step into the future:
   dealii::Vector<double>       old_solution;
   dealii::Vector<double>       current_solution;
   dealii::Vector<double>       predictor;
   
   dealii::Vector<double>       right_hand_side;
   
   // This final set of member variables
   // (except for the object holding all
   // run-time parameters at the very
   // bottom and a screen output stream
   // that only prints something if
   // verbose output has been requested)
   // deals with the inteface we have in
   // this program to the Trilinos library
   // that provides us with linear
   // solvers. Similarly to including
   // PETSc matrices in step-17,
   // step-18, and step-19, all we
   // need to do is to create a Trilinos
   // sparse matrix instead of the
   // standard deal.II class. The system
   // matrix is used for the Jacobian in
   // each Newton step. Since we do not
   // intend to run this program in
   // parallel (which wouldn't be too hard
   // with Trilinos data structures,
   // though), we don't have to think
   // about anything else like
   // distributing the degrees of freedom.
   dealii::TrilinosWrappers::SparseMatrix system_matrix;
   
   Parameters::AllParameters<dim>  parameters;
   dealii::ConditionalOStream      verbose_cout;

   // Call the appropriate numerical flux function
   template <typename InputVector>
   void numerical_normal_flux (const dealii::Point<dim>  &normal,
                               const InputVector         &Wplus,
                               const InputVector         &Wminus,
                               Sacado::Fad::DFad<double> (&normal_flux)[EulerEquations<dim>::n_components])
   {
      switch(parameters.flux_type)
      {
         case Parameters::Flux::lxf:
            EulerEquations<dim>::lxf_flux (normal,
                                           Wplus,
                                           Wminus,
                                           normal_flux);
            break;

         case Parameters::Flux::sw:
            EulerEquations<dim>::steger_warming_flux (normal,
                                                      Wplus,
                                                      Wminus,
                                                      normal_flux);
            break;

         case Parameters::Flux::kfvs:
            EulerEquations<dim>::kfvs_flux (normal,
                                            Wplus,
                                            Wminus,
                                            normal_flux);
            break;

	      default:
            Assert (false, dealii::ExcNotImplemented());
      }
   }

};
#endif
