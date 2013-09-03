#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <deal.II/base/logstream.h>

using namespace dealii;

double alpha, Omega, W, intA;

class InfinityValues : public Function<2>
{
   public:
        InfinityValues () : Function<2>() {}
          virtual double value (const Point<2>   &p,
                                const unsigned int  component = 0) const;
};

double InfinityValues::value (const Point<2> &p,
                              const unsigned int /*component*/) const
{
   double x = p(0);
   double y = p(1);
   return -0.5 * W * x * x + 0.25 * Omega * x * x * intA / std::pow(x*x + y*y, 1.5);
}

class Norbury
{
   public:
      Norbury ();
      void run ();

   private:
      void read_grid ();
      void setup_system ();
      void assemble_system ();
      void solve ();
      void output_results ();

      Triangulation<2>   triangulation;
      FE_Q<2>            fe;
      FE_DGQ<2>          fecell;
      DoFHandler<2>      dof_handler;
      DoFHandler<2>      dof_handler_cell;
      SparsityPattern      sparsity_pattern;
      SparseMatrix<double> system_matrix;
      Vector<double>       psi;
      Vector<double>       charfun;
      Vector<double>       system_rhs;
};

Norbury::Norbury ()
   :
      fe (1),
      fecell (0),
      dof_handler (triangulation),
      dof_handler_cell (triangulation)
{
}

void Norbury::read_grid ()
{
   GridIn<2> grid_in;
   grid_in.attach_triangulation(triangulation);
   std::ifstream input_file("norbury.msh");
   grid_in.read_msh(input_file);
}

void Norbury::setup_system ()
{
   dof_handler.distribute_dofs (fe);
   dof_handler_cell.distribute_dofs (fecell);
   std::cout << "   Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << std::endl;
   CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
   sparsity_pattern.copy_from(c_sparsity);
   system_matrix.reinit (sparsity_pattern);
   psi.reinit (dof_handler.n_dofs());
   system_rhs.reinit (dof_handler.n_dofs());

   charfun.reinit (dof_handler_cell.n_dofs());
   QGauss<2>  quadrature_formula(1);
   FEValues<2> fe_values (fecell, quadrature_formula,
                          update_JxW_values);
   std::vector<types::global_dof_index> local_dof_indices (1);
   
   typename DoFHandler<2>::active_cell_iterator
      cell = dof_handler_cell.begin_active(),
      endc = dof_handler_cell.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      if(cell->material_id() == 2)
         charfun(local_dof_indices[0]) = 1;
      else
         charfun(local_dof_indices[0]) = 0;
   }
}

void Norbury::assemble_system ()
{
   QGauss<2>  quadrature_formula(2);
   FEValues<2> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

   intA = 0;

   typename DoFHandler<2>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         const double x = fe_values.quadrature_point (q_point)(0);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                   fe_values.shape_grad (j, q_point) *
                                   fe_values.JxW (q_point)) / x;
            if(cell->material_id() == 2)
            {
               cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                              Omega * x *
                              fe_values.JxW (q_point));
               intA += x * x * x * fe_values.JxW (q_point);
            }
          }
      }


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
   }

   // We have only half the domain
   intA *= 2.0;
   std::cout << "intA = " << intA << std::endl;

  std::map<types::global_dof_index,double> boundary_values;
  // on x=0, we have psi=0
  VectorTools::interpolate_boundary_values (dof_handler,
                                            4,
                                            ZeroFunction<2>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      psi,
                                      system_rhs);

  // farfield
  VectorTools::interpolate_boundary_values (dof_handler,
                                            5,
                                            InfinityValues(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      psi,
                                      system_rhs);
}

void Norbury::solve ()
{
   SolverControl   solver_control (2000, 1e-12);
   SolverCG<>      solver (solver_control);
   PreconditionSSOR<> preconditioner;
   preconditioner.initialize(system_matrix, 1.2);
   solver.solve (system_matrix, psi, system_rhs,
                preconditioner);
   std::cout << "   " << solver_control.last_step()
             << " CG iterations needed to obtain convergence."
             << std::endl;
}

void Norbury::output_results ()
{
   DataOut<2> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (psi, "psi");
   data_out.build_patches ();
   std::ofstream output ("psi.vtk");
   data_out.write_vtk (output);

   DataOut<2> data_out_cell;
   data_out_cell.attach_dof_handler (dof_handler_cell);
   data_out_cell.add_data_vector (charfun, "charfun");
   data_out_cell.build_patches ();
   std::ofstream output_cell ("charfun.vtk");
   data_out_cell.write_vtk (output_cell);
}

void Norbury::run ()
{
   read_grid();
   setup_system ();
   assemble_system ();
   solve ();
   output_results ();
}

int main ()
{
   deallog.depth_console (0);
   W = 0.6586;
   alpha = 0.4;
   Omega = 1.0/alpha/alpha;
   Norbury norbury;
   norbury.run ();
   return 0;
}
