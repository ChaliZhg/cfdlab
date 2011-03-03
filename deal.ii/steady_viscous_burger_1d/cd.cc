#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

template <int dim>
class ConvDiffProblem
{
public:
    ConvDiffProblem ();
    void run ();

private:
    void make_grid_and_dofs ();
    void assemble_system_matrix ();
    void assemble_system_rhs ();
    void solve ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    double               dt;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};

// Source term
template <int dim>
class SourceTerm : public Function<dim>
{
public:
    SourceTerm () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};

// Source term value
template <int dim>
double SourceTerm<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
   double a   = 10.0;
   double x   = p[0];
   double u   = 10.0 * x * (1.0 - x) * sin(a * x); // Exact solution
   double ux  = 10.0 * a * (1.0-x) * x * cos(a * x) + 
                10.0 * (1.0 - x) * sin(a * x) -
                10.0 * x * sin(a * x);
   double uxx = -20.0 * (a * x * cos(a * x) + sin(a * x)) +
                 10.0 * (1.0-x) * (2.0 * a * cos(a * x) - a * a * x * sin(a * x));

   return u * ux - uxx;
}

// Boundary condition
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
    BoundaryValues () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


// Boundary condition value
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
   double a = 10.0;
   return 10.0 * p[0] * (1.0 - p[0]) * sin(a * p[0]);
}

// Constructor
template <int dim>
ConvDiffProblem<dim>::ConvDiffProblem () :
    fe (1),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
}

template <int dim>
void ConvDiffProblem<dim>::make_grid_and_dofs ()
{
    GridGenerator::hyper_cube (triangulation, 0, 1);
    triangulation.refine_global (4);

    std::cout << "   Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "   Total number of cells: "
              << triangulation.n_cells()
              << std::endl;

    dof_handler.distribute_dofs (fe);

    std::cout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
}

// Assemble mass matrix
template <int dim>
void ConvDiffProblem<dim>::assemble_system_matrix ()
{
    std::cout << "Constructing system matrix ...\n";

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
        cell_matrix = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (
                                         (fe_values.shape_value (i, q_point) *
                                          fe_values.shape_value (j, q_point) ) / dt
                                        +(fe_values.shape_grad (i, q_point) *
                                          fe_values.shape_grad (j, q_point) )
                                        ) * fe_values.JxW (q_point);
            }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
                system_matrix.add (local_dof_indices[i],
                                   local_dof_indices[j],
                                   cell_matrix(i,j));
    }

}
// Assemble system rhs
template <int dim>
void ConvDiffProblem<dim>::assemble_system_rhs ()
{
    QGauss<dim>  quadrature_formula(2);

    const SourceTerm<dim> source_term;

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  solution_values (n_q_points);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    // Initialize to zero since this function is called repeatedly
    system_rhs = 0.0;

    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
        cell_rhs = 0;

        fe_values.get_function_values (solution, solution_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            double flux = 0.5 * solution_values[q_point] * solution_values[q_point];
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                               solution_values[q_point]/dt *
                               fe_values.JxW (q_point));

                cell_rhs(i) += (fe_values.shape_grad (i, q_point)[0] *
                                flux *
                                fe_values.JxW (q_point));

                cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                                source_term.value (fe_values.quadrature_point (q_point)) *
                                fe_values.JxW (q_point));
            }
        }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

	// left boundary condition
    std::map<unsigned int,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
            0,
            BoundaryValues<dim>(),
            boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
	// right boundary condition
    VectorTools::interpolate_boundary_values (dof_handler,
            1,
            BoundaryValues<dim>(),
            boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
}

// Update solution by one time step
template <int dim>
void ConvDiffProblem<dim>::solve ()
{
    SolverControl           solver_control (1000, 1e-4);
    SolverCG<>              cg (solver_control);
    cg.solve (system_matrix, solution, system_rhs,
              PreconditionIdentity());

    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence."
              << std::endl;
}

// Save solution to file
template <int dim>
void ConvDiffProblem<dim>::output_results () const
{
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");

    data_out.build_patches ();

    std::ofstream output ("solution.gpl");
    data_out.write_gnuplot (output);
}

// Start solving the problem
template <int dim>
void ConvDiffProblem<dim>::run ()
{
    std::cout << "Solving convection-diffusion problem ...\n";

    make_grid_and_dofs();

    dt = 1.0e-2;
    assemble_system_matrix ();


    unsigned int iter = 0;
    solution = 0.0;
    while (iter < 20)
    {
      assemble_system_rhs ();
      solve ();
      ++iter;
    }
    output_results ();
}

// Main function
int main ()
{
    deallog.depth_console (0);
    {
        ConvDiffProblem<1> cd_problem;
        cd_problem.run ();
    }

    return 0;
}

