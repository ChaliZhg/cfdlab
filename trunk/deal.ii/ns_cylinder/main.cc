/*
 Solution of incompressible Navier-Stokes equations using Taylor-Hood
 
 BDF1 in first time step; BDF2 from second time step onwards
 Convective term is linearized by extrapolation
 Second order accurate in time
 
 Author: Praveen. C, 
         TIFR-CAM, Bangalore
         http://praveen.tifrbng.res.in
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>(dim+1) {}
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
   virtual void vector_value (const Point<dim> &p,
                              Vector<double>   &value) const;
};

template <int dim>
double
InitialCondition<dim>::value (const Point<dim>  &p,
                            const unsigned int component) const
{
   Assert (component < this->n_components,
           ExcIndexRange (component, 0, this->n_components));
   if (component == 0)
      return (1.5/std::pow(0.205,2))*(p[1]+0.2)*(0.21-p[1]);
   else
      return 0.0;
}

template <int dim>
void
InitialCondition<dim>::vector_value (const Point<dim> &p,
                                   Vector<double>   &values) const
{
   for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = InitialCondition<dim>::value (p, c);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
class NS
{
public:
   NS (unsigned int degree);
   ~NS() {};
   void run ();
   
private:
   void run_steady ();
   void run_unsteady ();
   void make_grid_dofs ();
   void assemble_mass_matrix ();
   void assemble_matrix (unsigned int order);
   void assemble_matrix_and_rhs (unsigned int order);
   void solve ();
   void compute_vorticity ();
   void output_results() const;
   
   unsigned int               degree;
   FESystem<dim>              fe;
   FE_Q<dim>                  fe_scalar;
   Triangulation<dim>         triangulation;
   DoFHandler<dim>            dof_handler;
   DoFHandler<dim>            dof_handler_scalar;
   MappingQ<dim>              mapping;
   
   ConstraintMatrix           constraints;
   BlockSparsityPattern       sparsity_pattern;
   BlockSparseMatrix<double>  system_matrix_constant;
   BlockSparseMatrix<double>  system_matrix;
   BlockVector<double>        solution0, solution1, solution2;
   BlockVector<double>        system_rhs;
   
   SparsityPattern            sparsity_pattern_scalar;
   SparseMatrix<double>       mass_matrix;
   Vector<double>             vorticity;
   SparseDirectUMFPACK        vorticity_solver;
   
   // Parameters
   double                     dt, Uref, Lref, Re, viscosity, final_time;
};

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
NS<dim>::NS (unsigned int degree)
:
   degree (degree),
   fe( FE_Q<dim>(QGaussLobatto<1>(degree+2)), dim,
       FE_Q<dim>(QGaussLobatto<1>(degree+1)),   1),
   fe_scalar (FE_Q<dim>(QGaussLobatto<1>(degree+2))),
   dof_handler (triangulation),
   dof_handler_scalar (triangulation),
   mapping (degree+1)
{
   dt = 0.01;
   Re = 100.0;
   Uref = 1.0;
   Lref = 0.1;
   viscosity = Uref*Lref/Re;
   final_time = 40.0;
   std::string grid_file = "karman.msh";
   
   std::cout << "Reading grid from " << grid_file << std::endl;
   GridIn<dim> grid_in;
   grid_in.attach_triangulation (triangulation);
   std::ifstream input_file (grid_file.c_str());
   grid_in.read_msh (input_file);
   
   double radius = Lref / 2.0;
   Point<dim> center (0.0, 0.0);
   static const HyperBallBoundary<dim> boundary_description (center, radius);
   triangulation.set_boundary (2, boundary_description);
   
   std::ofstream grid_output_file("grid.eps");
   GridOut grid_out;
   grid_out.write_eps (triangulation, grid_output_file);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::make_grid_dofs()
{
   dof_handler.distribute_dofs (fe);
   
   DoFRenumbering::Cuthill_McKee (dof_handler);
   std::vector<unsigned int> block_component (dim+1,0);
   block_component[dim] = 1;
   DoFRenumbering::component_wise (dof_handler, block_component);

   std::vector<types::global_dof_index> dofs_per_block (2);
   DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
   const unsigned int n_u = dofs_per_block[0],
                      n_p = dofs_per_block[1];
   std::cout << "   Number of active cells: "
             << triangulation.n_active_cells()
             << std::endl
             << "   Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << " (" << n_u << '+' << n_p << ')'
             << std::endl;
   
   {
      BlockCompressedSimpleSparsityPattern csp (2,2);
      csp.block(0,0).reinit (n_u, n_u);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,1).reinit (n_p, n_p);
      csp.collect_sizes();
      DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
      sparsity_pattern.copy_from (csp);
   }
   
   system_matrix_constant.reinit (sparsity_pattern);
   system_matrix.reinit (sparsity_pattern);
   
   solution0.reinit (2);
   solution0.block(0).reinit (n_u);
   solution0.block(1).reinit (n_p);
   solution0.collect_sizes ();

   solution1.reinit (2);
   solution1.block(0).reinit (n_u);
   solution1.block(1).reinit (n_p);
   solution1.collect_sizes ();

   solution2.reinit (2);
   solution2.block(0).reinit (n_u);
   solution2.block(1).reinit (n_p);
   solution2.collect_sizes ();

   system_rhs.reinit (2);
   system_rhs.block(0).reinit (n_u);
   system_rhs.block(1).reinit (n_p);
   system_rhs.collect_sizes ();
   
   // These are needed for computing vorticity
   dof_handler_scalar.distribute_dofs (fe_scalar);
   DoFRenumbering::Cuthill_McKee (dof_handler_scalar);

   {
      CompressedSparsityPattern csp (dof_handler_scalar.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler_scalar, csp);
      sparsity_pattern_scalar.copy_from (csp);
   }
   
   mass_matrix.reinit (sparsity_pattern_scalar);
   vorticity.reinit (dof_handler_scalar.n_dofs());
   std::cout << "   Number of vorticity dofs: "
             << dof_handler_scalar.n_dofs()
             << std::endl;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_mass_matrix ()
{
   mass_matrix = 0;
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe_scalar, quadrature_formula,
                            update_values    |
                            update_JxW_values);
   const unsigned int   dofs_per_cell   = fe_scalar.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler_scalar.begin_active(),
      endc = dof_handler_scalar.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      local_matrix = 0;
      
      for (unsigned int q=0; q<n_q_points; ++q)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<=i; ++j)
               local_matrix(i,j) +=   fe_values.shape_value(i,q)
                                    * fe_values.shape_value(j,q)
                                    * fe_values.JxW(q);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            mass_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
   }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_matrix (unsigned int order)
{
   double a2;
   if(order == 0)
      a2 = 0.0;
   else if(order == 1)
      a2 = 1.0;
   else if(order == 2)
      a2 = 1.5;
   else
      Assert(false, ExcMessage("Not implemented"));
   
   system_matrix_constant = 0;
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_values    |
                            update_gradients |
                            update_JxW_values);
   const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   const FEValuesExtractors::Scalar pressure (dim);
   
   std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
   std::vector<double>                  div_phi_u     (dofs_per_cell);
   std::vector<double>                  phi_p         (dofs_per_cell);
   std::vector<Tensor<1,dim> >          phi_u         (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      local_matrix = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
      {
         for (unsigned int k=0; k<dofs_per_cell; ++k)
         {
            symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
            div_phi_u[k]     = fe_values[velocities].divergence (k, q);
            phi_p[k]         = fe_values[pressure].value (k, q);
            phi_u[k]         = fe_values[velocities].value (k, q);
         }
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<=i; ++j)
            {
               local_matrix(i,j) += ((a2/dt) * phi_u[i] * phi_u[j]
                                     + 2.0 * viscosity * symgrad_phi_u[i] * symgrad_phi_u[j]
                                     - div_phi_u[i] * phi_p[j]
                                     - phi_p[i] * div_phi_u[j]
                                     )
                                     * fe_values.JxW(q);
            }
         }
      }
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix_constant.add (local_dof_indices[i],
                                        local_dof_indices[j],
                                        local_matrix(i,j));
   }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_matrix_and_rhs (unsigned int order)
{
   double a0, a1;
   if(order == 0)
   {
      a0 = 0.0;
      a1 = 0.0;
   }
   else if(order == 1)
   {
      a0 =-1.0;
      a1 = 0.0;
   }
   else if(order == 2)
   {
      a0 =  0.5;
      a1 = -2.0;
   }
   else
      Assert(false, ExcMessage("Not implemented"));
   
   // use solution2 for extrapolated velocity
   if(order == 0)
      ; // nothing to do; we use solution2
   else if(order == 1)
      solution2.block(0) = solution0.block(0);
   else
      // solution2 = 2 * solution1 - solution0
      solution2.block(0).equ(2.0, solution1.block(0), -1.0, solution0.block(0));
   
   system_matrix.copy_from(system_matrix_constant);
   system_rhs    = 0;
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_values    |
                            update_gradients |
                            update_JxW_values);
   const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       local_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   const FEValuesExtractors::Scalar pressure (dim);
   
   std::vector<Tensor<2,dim> >  grad_phi_u (dofs_per_cell);
   std::vector<Tensor<1,dim> >  phi_u         (dofs_per_cell);
   std::vector<Tensor<1,dim> >  velocity(n_q_points, Tensor<1,dim>());
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      fe_values[velocities].get_function_values (solution2, velocity);
      local_matrix = 0;
      local_rhs    = 0;
      cell->get_dof_indices (local_dof_indices);
      
      for (unsigned int q=0; q<n_q_points; ++q)
      {
         for (unsigned int k=0; k<dofs_per_cell; ++k)
         {
            grad_phi_u[k] = fe_values[velocities].gradient (k, q);
            phi_u[k]      = fe_values[velocities].value (k, q);
         }
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
               local_matrix(i,j) += (grad_phi_u[j] * velocity[q]) * phi_u[i]
                                     * fe_values.JxW(q);
               local_rhs(i) += (1.0/dt) *
                               (-a1 * solution1(local_dof_indices[j]) - a0 * solution0(local_dof_indices[j]))
                               * phi_u[i] * phi_u[j] * fe_values.JxW(q);
            }
         }
      }
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         system_rhs(local_dof_indices[i]) += local_rhs(i);
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               local_matrix(i,j));
      }
   }
   
   // Apply boundary conditions
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             1,
                                             InitialCondition<dim>(),
                                             boundary_values,
                                             fe.component_mask(velocities));
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             2,
                                             ZeroFunction<dim>(dim+1),
                                             boundary_values,
                                             fe.component_mask(velocities));
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             3,
                                             ZeroFunction<dim>(dim+1),
                                             boundary_values,
                                             fe.component_mask(velocities));
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution2,
                                       system_rhs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::solve()
{
   SparseDirectUMFPACK  solver;
   solver.initialize (system_matrix);
   solver.vmult (solution2, system_rhs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::compute_vorticity ()
{
   static unsigned int status = 0;
   
   if(status == 0)
   {
      assemble_mass_matrix ();
      vorticity_solver.initialize (mass_matrix);
      status = 1;
   }
   
   Vector<double> vorticity_rhs (dof_handler_scalar.n_dofs());
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_gradients);
   FEValues<dim> fe_values_vorticity (mapping, fe_scalar, quadrature_formula,
                            update_values | update_JxW_values);
   const unsigned int   dofs_per_cell   = fe_scalar.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   Vector<double>       local_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   std::vector<typename FEValuesViews::Vector<dim>::curl_type> vorticity_values (n_q_points);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell_vorticity = dof_handler_scalar.begin_active();
   for (; cell!=endc; ++cell, ++cell_vorticity)
   {
      fe_values.reinit (cell);
      fe_values_vorticity.reinit (cell_vorticity);

      fe_values[velocities].get_function_curls (solution2, vorticity_values);
      
      local_rhs    = 0;
      
      for (unsigned int q=0; q<n_q_points; ++q)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            local_rhs(i) +=   vorticity_values[q][0]
                            * fe_values_vorticity.shape_value(i,q)
                            * fe_values_vorticity.JxW(q);
      
      cell_vorticity->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         vorticity_rhs(local_dof_indices[i]) += local_rhs(i);
   }
   
   vorticity_solver.vmult (vorticity, vorticity_rhs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void
NS<dim>::output_results ()  const
{
   static unsigned int cycle = 0;

   std::vector<std::string> solution_names (dim, "velocity");
   solution_names.push_back ("pressure");
   
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
   data_component_interpretation
   (dim, DataComponentInterpretation::component_is_part_of_vector);
   data_component_interpretation
   .push_back (DataComponentInterpretation::component_is_scalar);
   
   DataOut<dim> data_out;
   //data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (dof_handler, solution2, solution_names,
                             /*DataOut<dim>::type_dof_data,*/
                             data_component_interpretation);
   data_out.add_data_vector (dof_handler_scalar, vorticity, "vorticity"/*,
                             DataOut<dim>::type_dof_data*/);

   data_out.build_patches (mapping, degree+1);
   
   std::ostringstream filename;
   filename << "solution-"
            << Utilities::int_to_string (cycle, 3)
            << ".vtk";
   std::ofstream output (filename.str().c_str());
   data_out.write_vtk (output);
   
   ++cycle;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run_steady ()
{
   unsigned int order = 0;
   assemble_matrix (order);
   
   // Set initial condition
   std::cout << "Setting initial condition ..." << std::endl;
   VectorTools::interpolate(mapping, dof_handler,
                            InitialCondition<dim>(), solution0);
   solution1 = solution0;
   solution2 = solution0;
   compute_vorticity ();
   output_results ();
   
   unsigned int iter = 0;
   
   while (iter < 10)
   {
      // Assemble matrix and rhs
      assemble_matrix_and_rhs (order);
      
      // solve
      solve ();
      
      ++iter;
      std::cout << iter << std::endl;
      
      compute_vorticity ();
      output_results ();
   }
   
   // save solution to file
   std::ofstream output_file("steady.dat");
   solution2.block_write (output_file);
   
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run_unsteady ()
{
   unsigned int order = 1;
   assemble_matrix (order);
   
   
   // Set initial condition
   std::cout << "Setting initial condition ..." << std::endl;
   VectorTools::interpolate(mapping, dof_handler,
                            InitialCondition<dim>(), solution0);
   solution1 = solution0;
   solution2 = solution0;
   compute_vorticity ();
   output_results ();
   
   double time = 0;
   unsigned int iter = 0;
   
   while (time < final_time)
   {
      // Assemble matrix and rhs
      assemble_matrix_and_rhs (order);
      
      // solve
      solve ();
      
      time += dt;
      ++iter;
      std::cout << iter << "  " << time << "  " << std::endl;
      
      if(iter == 1)
      {
         solution1 = solution2;

         order = 2;
         assemble_matrix (order);
      }
      else
      {
         solution0 = solution1;
         solution1 = solution2;
      }
      
      if(iter%10 == 0)
      {
         compute_vorticity ();
         output_results ();
      }
   }
   
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run ()
{
   make_grid_dofs ();

   //run_steady ();
   run_unsteady ();
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   try
   {
      deallog.depth_console (0);

      unsigned int pressure_degree = 2;
      NS<2> ns_problem (pressure_degree);
      ns_problem.run();
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
   }
   return 0;
}