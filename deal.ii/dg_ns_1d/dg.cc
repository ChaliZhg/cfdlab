#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_dgq.h>
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

// Number of variables: mass, momentum and energy
const unsigned int n_var = 3;
const double gas_gamma = 1.4;

// Function declaration
std::vector<double> ns_flux (const double density,
                             const double momentum,
                             const double energy);

// Main class of the problem
template <int dim>
class NSProblem
{
public:
    NSProblem (unsigned int degree);
    void run ();

private:
    void make_grid_and_dofs ();
    void initialize ();
    void assemble_mass_matrix ();
    void assemble_rhs ();
    void update ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> inv_mass_matrix;

    double               dt;
    Vector<double>       density;
    Vector<double>       momentum;
    Vector<double>       energy;
    Vector<double>       rhs_density;
    Vector<double>       rhs_momentum;
    Vector<double>       rhs_energy;
};

// Constructor
template <int dim>
NSProblem<dim>::NSProblem (unsigned int degree) :
    fe (degree),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
}

// Make grid and allocate memory for solution variables
template <int dim>
void NSProblem<dim>::make_grid_and_dofs ()
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
   
    inv_mass_matrix.reinit (sparsity_pattern);
   
    // Solution variables
    density.reinit (dof_handler.n_dofs());
    rhs_density.reinit (dof_handler.n_dofs());
   
    momentum.reinit (dof_handler.n_dofs());
    rhs_momentum.reinit (dof_handler.n_dofs());
   
    energy.reinit (dof_handler.n_dofs());
    rhs_energy.reinit (dof_handler.n_dofs());   
}

// Set initial conditions
template <int dim>
void NSProblem<dim>::initialize ()
{
}

// Assemble mass matrix for each cell
// Invert it and store
template <int dim>
void NSProblem<dim>::assemble_mass_matrix ()
{
    std::cout << "Constructing mass matrix ...\n";

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   inv_cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
    // Cell iterator
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
        cell_matrix = 0.0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += fe_values.shape_value (i, q_point) *
                                        fe_values.shape_value (j, q_point) *
                                        fe_values.JxW (q_point);
            }
       
        // Invert cell_matrix
        inv_cell_matrix.invert(cell_matrix);
       
        // Store the inverse
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
           for (unsigned int j=0; j<dofs_per_cell; ++j)
              inv_mass_matrix.set (local_dof_indices[i],
                                   local_dof_indices[j],
                                   inv_cell_matrix(i,j));          
    }

}

// Flux for NS equation
std::vector<double> ns_flux (const double density,
                             const double momentum,
                             const double energy)
{
   std::vector<double> flux(3);
   
   double velocity = momentum / density;
   double pressure = (gas_gamma - 1.0) * (energy - 0.5 * momentum / density);
   flux[0] = momentum;
   flux[1] = pressure + momentum * velocity;
   flux[2] = (energy + pressure) * velocity;
   
   return flux;
}

// Assemble system rhs
template <int dim>
void NSProblem<dim>::assemble_rhs ()
{
    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   |
                             update_quadrature_points | 
                             update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  density_values (n_q_points);
    std::vector<double>  momentum_values (n_q_points);
    std::vector<double>  energy_values (n_q_points);

    Vector<double>       cell_rhs_density  (dofs_per_cell);
    Vector<double>       cell_rhs_momentum (dofs_per_cell);
    Vector<double>       cell_rhs_energy   (dofs_per_cell);


    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
       
        cell_rhs_density  = 0.0;
        cell_rhs_momentum = 0.0;
        cell_rhs_energy   = 0.0;

        // Compute conserved variables at quadrature points
        fe_values.get_function_values (density,  density_values);
        fe_values.get_function_values (momentum, momentum_values);
        fe_values.get_function_values (energy,   energy_values);

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            std::vector<double> flux = ns_flux(density_values[q_point], 
                                               momentum_values[q_point], 
                                               energy_values[q_point]);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                cell_rhs_density(i) += (fe_values.shape_grad (i, q_point)[0] *
                                        flux[0] *
                                        fe_values.JxW (q_point));
                cell_rhs_momentum(i)+= (fe_values.shape_grad (i, q_point)[0] *
                                        flux[1] *
                                        fe_values.JxW (q_point));
                cell_rhs_energy(i)  += (fe_values.shape_grad (i, q_point)[0] *
                                        flux[2] *
                                        fe_values.JxW (q_point));
            }
        }

        // Compute flux at cell boundaries TBD
       

       
        // Multiply by inverse mass matrix and add to rhs
        cell->get_dof_indices (local_dof_indices);
        unsigned int ig, jg;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            ig = local_dof_indices[i];
           
            rhs_density (ig) = 0.0;
            rhs_momentum(ig) = 0.0;
            rhs_energy  (ig) = 0.0;
           
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
               jg = local_dof_indices[j];
               rhs_density(ig)  += inv_mass_matrix(ig,jg) * cell_rhs_density(j);
               rhs_momentum(ig) += inv_mass_matrix(ig,jg) * cell_rhs_momentum(j);
               rhs_energy(ig)   += inv_mass_matrix(ig,jg) * cell_rhs_energy(j);
            }

        }
    }

}

// Update solution by one time step
// Forward Euler; TBD implement 3-stage RK
template <int dim>
void NSProblem<dim>::update ()
{
   // Multiply by time step
   rhs_density.scale(dt);
   rhs_momentum.scale(dt);
   rhs_energy.scale(dt);
   
   // Update conserved variables
   density  += rhs_density;
   momentum += rhs_momentum;
   energy   += rhs_energy;
}

// Save solution to file
template <int dim>
void NSProblem<dim>::output_results () const
{
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (density, "density");

    data_out.build_patches ();

    std::ofstream output ("solution.gpl");
    data_out.write_gnuplot (output);
}

// Start solving the problem
template <int dim>
void NSProblem<dim>::run ()
{
    std::cout << "Solving 1-D NS problem ...\n";

    make_grid_and_dofs();

    dt = 1.0e-2;
    assemble_mass_matrix ();
    initialize ();

    unsigned int iter = 0;
    while (iter < 20)
    {
      assemble_rhs ();
      update ();
      ++iter;
    }
    output_results ();
}

// Main function
int main ()
{
    deallog.depth_console (0);
    {
        NSProblem<1> ns_problem(0);
        ns_problem.run ();
    }

    return 0;
}

