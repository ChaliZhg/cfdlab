/* $Id: step-12.cc 22886 2010-11-29 23:48:38Z bangerth $ */
/* Author: Guido Kanschat, Texas A&M University, 2009 */

/*    $Id: step-12.cc 22886 2010-11-29 23:48:38Z bangerth $       */
/*                                                                */
/*    Copyright (C) 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

// Modifications by Praveen. C, http://math.tifrbng.res.in/~praveen
// Explicit time-stepping Runge-Kutta DG method
// Mass matrix on each cell is computed, inverted and the inverse 
// is stored. Then in each time iteration, we need to compute right
// hand side and multipy by inverse mass mass matrix. After that
// solution is advanced to new time level by an RK scheme.

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe_values.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/data_out.h>
#include <fe/mapping_q1.h>
		
#include <fe/fe_dgq.h>
				 
#include <numerics/derivative_approximation.h>

#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_info.h>
#include <numerics/mesh_worker_assembler.h>
#include <numerics/mesh_worker_loop.h>

#include <iostream>
#include <fstream>

using namespace dealii;

//------------------------------------------------------------------------------
// Boundary condition function class
//------------------------------------------------------------------------------
template <int dim>
class BoundaryValues: public Function<dim>
{
  public:
    BoundaryValues () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
			                    std::vector<double> &values,
			                    const unsigned int component=0) const;
};

// Computes boundary condition value at a list of boundary points
template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
				       std::vector<double> &values,
				       const unsigned int) const
{
   Assert(values.size()==points.size(),
          ExcDimensionMismatch(values.size(),points.size()));
   
   for (unsigned int i=0; i<values.size(); ++i)
   {
      if (points[i](0)<0.5)
         values[i]=1.;
      else
         values[i]=0.;
   }
}

//------------------------------------------------------------------------------
// Class for integrating rhs using MeshWorker
//------------------------------------------------------------------------------
template <int dim>
class RHSIntegrator
{
   public:
      RHSIntegrator (const DoFHandler<dim>& dof_handler)
         : dof_info (dof_handler) {};

      MeshWorker::IntegrationInfoBox<dim> info_box;
      MeshWorker::DoFInfo<dim> dof_info;
      MeshWorker::Assembler::ResidualSimple< Vector<double> >
         assembler;
};

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class Step12
{
   public:
      Step12 ();
      void run ();
      
   private:
      void setup_system ();
      void assemble_mass_matrix ();
      void setup_mesh_worker (RHSIntegrator<dim>&);
      void assemble_rhs (RHSIntegrator<dim>&);
      void solve ();
      void refine_grid ();
      void output_results (const unsigned int cycle) const;
      
      Triangulation<dim>   triangulation;
      const MappingQ1<dim> mapping;
      
      FE_DGQ<dim>          fe;
      DoFHandler<dim>      dof_handler;
      
      std::vector< FullMatrix<double> > inv_mass_matrix;
      
      Vector<double>       solution;
      Vector<double>       solution_old;
      Vector<double>       right_hand_side;
   
      typedef MeshWorker::DoFInfo<dim> DoFInfo;
      typedef MeshWorker::IntegrationInfo<dim> CellInfo;
      
      static void integrate_cell_term (DoFInfo& dinfo, CellInfo& info);
      static void integrate_boundary_term (DoFInfo& dinfo, CellInfo& info);
      static void integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                       CellInfo& info1, CellInfo& info2);
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
Step12<dim>::Step12 ()
      :
      mapping (),
      fe (1),
      dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
// Make dofs and allocate memory
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_system ()
{
   std::cout << "Allocating memory ...\n";

   dof_handler.distribute_dofs (fe);
   
   inv_mass_matrix.resize (triangulation.n_cells(), 
                           FullMatrix<double>(fe.dofs_per_cell,
                                              fe.dofs_per_cell));
   solution.reinit (dof_handler.n_dofs());
   solution_old.reinit (dof_handler.n_dofs());
   right_hand_side.reinit (dof_handler.n_dofs());
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_mass_matrix ()
{
   std::cout << "Constructing mass matrix ...\n";
   
   QGauss<dim>  quadrature_formula(fe.degree+1);
   
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values | update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
      
   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      fe_values.reinit (cell);
      cell_matrix = 0.0;
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
               cell_matrix(i,j) += fe_values.shape_value (i, q_point) *
                                   fe_values.shape_value (j, q_point) *
                                   fe_values.JxW (q_point);
      
      // Invert cell_matrix
      inv_mass_matrix[c].invert (cell_matrix);      
   }
   
}

//------------------------------------------------------------------------------
// Create mesh worker for integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_mesh_worker (RHSIntegrator<dim>& rhs_integrator)
{   
   std::cout << "Setting up mesh worker ...\n";

   MeshWorker::IntegrationInfoBox<dim>& info_box = rhs_integrator.info_box;
   MeshWorker::DoFInfo<dim>& dof_info = rhs_integrator.dof_info;
   MeshWorker::Assembler::ResidualSimple< Vector<double> >&
      assembler = rhs_integrator.assembler;

   const unsigned int n_gauss_points = dof_handler.get_fe().degree+1;
   info_box.initialize_gauss_quadrature(n_gauss_points,
                                        n_gauss_points,
                                        n_gauss_points);

   // Add solution vector to info_box
   NamedData< Vector<double>* > solution_data;
   solution_data.add (&solution, "solution");
   info_box.cell_selector.add     ("solution", true, false, false);
   info_box.boundary_selector.add ("solution", true, false, false);
   info_box.face_selector.add     ("solution", true, false, false);
   
   info_box.initialize_update_flags ();
   info_box.add_update_flags_all      (update_quadrature_points);
   info_box.add_update_flags_cell     (update_gradients);
   info_box.add_update_flags_boundary (update_values);
   info_box.add_update_flags_face     (update_values);

   info_box.initialize (fe, mapping, solution_data);
   
   // Attach rhs vector to assembler
   NamedData< Vector<double>* > rhs;
   rhs.add (&right_hand_side, "RHS");
   assembler.initialize (rhs);
}

//------------------------------------------------------------------------------
// Assemble rhs of the problem
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_rhs (RHSIntegrator<dim>& rhs_integrator)
{
   right_hand_side = 0.0;

   MeshWorker::integration_loop<dim, dim>
      (dof_handler.begin_active(), dof_handler.end(),
       rhs_integrator.dof_info, 
       rhs_integrator.info_box,
       &Step12<dim>::integrate_cell_term,
       &Step12<dim>::integrate_boundary_term,
       &Step12<dim>::integrate_face_term,
       rhs_integrator.assembler, true);

   // Multiply by inverse mass matrix
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   Vector<double> rhs (dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      cell->get_dof_indices (local_dof_indices);

      rhs = 0.0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            rhs(i) += inv_mass_matrix[c](i,j) * 
                      right_hand_side(local_dof_indices[j]);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
         right_hand_side(local_dof_indices[i]) = rhs(i);
   }
}

//------------------------------------------------------------------------------
// Compute cell integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_cell_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      beta(0) = -fe_v.quadrature_point(point)(1);
      beta(1) =  fe_v.quadrature_point(point)(0);
      beta /= beta.norm();
      
      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
         local_vector(i) += beta * fe_v.shape_grad(i,point) *
                            sol[point] * JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute boundary integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_boundary_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];

   Vector<double>& local_vector = dinfo.vector(0).block(0);
   
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Point<dim> >& normals = fe_v.get_normal_vectors ();
   
   std::vector<double> g(fe_v.n_quadrature_points);
   
   static BoundaryValues<dim> boundary_function;
   boundary_function.value_list (fe_v.get_quadrature_points(), g);
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      beta(0) = -fe_v.quadrature_point(point)(1);
      beta(1) = fe_v.quadrature_point(point)(0);
      beta /= beta.norm();
      
      const double beta_n=beta * normals[point];
      if (beta_n>0)
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector(i) -= beta_n * 
                               sol[point] *
                               fe_v.shape_value(i,point) *
                               JxW[point];
      else
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector(i) -= beta_n *
                               g[point] *
                               fe_v.shape_value(i,point) *
                               JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute integral over internal faces
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
				                           CellInfo& info1, CellInfo& info2)
{
   const FEValuesBase<dim>& fe_v          = info1.fe_values();
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();
   
   const std::vector<double>& sol1 = info1.values[0][0];
   const std::vector<double>& sol2 = info2.values[0][0];
   
   Vector<double>& local_vector1 = dinfo1.vector(0).block(0);
   Vector<double>& local_vector2 = dinfo2.vector(0).block(0);
   
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Point<dim> >& normals = fe_v.get_normal_vectors ();
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      beta(0) = -fe_v.quadrature_point(point)(1);
      beta(1) =  fe_v.quadrature_point(point)(0);
      beta   /= beta.norm();
      
      const double beta_n=beta * normals[point];
      if (beta_n>0)
      {
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= beta_n *
                                sol1[point] *
                                fe_v.shape_value(i,point) *
                                JxW[point];
         
         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += beta_n *
                                sol1[point] *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
      }
      else
      {
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= beta_n *
                                sol2[point] *
                                fe_v.shape_value(i,point) *
                                JxW[point];
         
         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += beta_n *
                                sol2[point] *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
      }
   }
}

//------------------------------------------------------------------------------
// Solve the problem to convergence by RK time integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::solve ()
{
   RHSIntegrator<dim> rhs_integrator (dof_handler);
   setup_mesh_worker (rhs_integrator);
   
   std::cout << "Solving by RK ...\n";

   double dt = 0.01;
   double residue0;
   double residue = 1.0e20;
   unsigned int iter = 0;
   while (iter < 1000 && residue > 1.0e-8)
   {
      assemble_rhs (rhs_integrator);

      for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
         solution(i) += dt * right_hand_side(i);

      residue = right_hand_side.l2_norm ();
      if(iter==0) residue0 = residue;
      residue /= residue0;
      
      ++iter;
   }
   
   std::cout << "Iterations=" << iter << ", Residue=" << residue << endl;

}

//------------------------------------------------------------------------------
// Refine grid
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::refine_grid ()
{
   Vector<float> gradient_indicator (triangulation.n_active_cells());
   
   DerivativeApproximation::approximate_gradient (mapping,
                                                  dof_handler,
                                                  solution,
                                                  gradient_indicator);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);
   
   GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                    gradient_indicator,
                                                    0.3, 0.1);
   
   triangulation.execute_coarsening_and_refinement ();
}

//------------------------------------------------------------------------------
// Save results to file
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::output_results (const unsigned int cycle) const
{
   // Write the grid in eps format.
   std::string filename = "grid-";
   filename += ('0' + cycle);
   Assert (cycle < 10, ExcInternalError());
   
   filename += ".eps";
   std::cout << "Writing grid to <" << filename << ">" << std::endl;
   std::ofstream eps_output (filename.c_str());
   
   GridOut grid_out;
   grid_out.write_eps (triangulation, eps_output);
   
   // Output of the solution in
   // gnuplot format.
   filename = "sol-";
   filename += ('0' + cycle);
   Assert (cycle < 10, ExcInternalError());
   
   filename += ".gnuplot";
   std::cout << "Writing solution to <" << filename << ">" << std::endl;
   std::ofstream gnuplot_output (filename.c_str());
   
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "u");
   
   data_out.build_patches ();
   
   data_out.write_gnuplot(gnuplot_output);
}

//------------------------------------------------------------------------------
// Actual computation starts from here
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::run ()
{
   for (unsigned int cycle=0; cycle<1; ++cycle)
   {
      std::cout << "Cycle " << cycle << std::endl;
      
      if (cycle == 0)
      {
         GridGenerator::hyper_cube (triangulation);
         triangulation.refine_global (5);
      }
      else
         refine_grid ();
      
      
      std::cout << "Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;
      
      setup_system ();
      
      std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;
      
      assemble_mass_matrix ();
      solve ();
      output_results (cycle);
   }
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main ()
{
   try
   {
      Step12<2> dgmethod;
      dgmethod.run ();
   }
   catch (std::exception &exc)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
		          << std::endl;
      std::cerr << "Exception on processing: " << std::endl
		          << exc.what() << std::endl
		          << "Aborting!" << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   };
   
   return 0;
}


