#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <numerics/matrix_tools.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <numerics/vector_tools.h>
#include <base/function.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <numerics/data_out.h>

#include <fstream>

using namespace dealii;

int main()
{
   Triangulation<2> triangulation;
   GridGenerator::hyper_cube (triangulation);
   triangulation.refine_global (4);

   const FE_Q<2> fe(1);

   DoFHandler<2> dof_handler (triangulation);
   dof_handler.distribute_dofs (fe);

   CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
   SparsityPattern sparsity_pattern;
   sparsity_pattern.copy_from(c_sparsity);

   SparseMatrix<double> system_matrix;
   system_matrix.reinit (sparsity_pattern);
   MatrixCreator::create_laplace_matrix (dof_handler, 
                                         QGauss<2>(2), 
                                         system_matrix);

   Vector<double> system_rhs;
   system_rhs.reinit (dof_handler.n_dofs());
   VectorTools::create_right_hand_side (dof_handler,
                                        QGauss<2>(2),
                                        ConstantFunction<2>(1.0),
                                        system_rhs);

   std::map<unsigned int,double> boundary_values;
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             ZeroFunction<2>(),
                                             boundary_values);

   Vector<double> solution;
   solution.reinit (dof_handler.n_dofs());

   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

   SolverControl solver_control (1000, 1.0e-12);
   SolverCG<>    solver (solver_control);
   solver.solve (system_matrix,
                 solution,
                 system_rhs,
                 PreconditionIdentity());

   DataOut<2> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches ();

   std::ofstream output ("solution.vtk");
   data_out.write_vtk (output);
}
