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

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

// Number of variables: mass, momentum and energy
const unsigned int n_var = 3;
const double gas_gamma = 1.4;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Initial condition
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>() {}
   
   virtual void vector_value (const Point<dim>   &p,
                              Vector<double>& values) const;
};

// Initial condition for density, velocity, pressure
template<int dim>
void InitialCondition<dim>::vector_value (const Point<dim>   &p,
                                          Vector<double>& values) const
{
   if(p[0] < 0.5)
   {
      values(0) = 1.0; // left density
      values(1) = 0.0; // left velocity
      values(2) = 1.0; // left pressure
   }
   else
   {
      values(0) = 0.125; // right density
      values(1) = 0.0; // right velocity
      values(2) = 0.1; // right pressure
   }

}

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
    void compute_face_flux ();
    void assemble_rhs ();
    void update (const unsigned int rk_stage);
    void output_results () const;
   
    double               xmin, xmax;
    unsigned int         n_cells;
    double               dt;
    double               dx;
    unsigned int         n_rk_stages;

    Triangulation<dim>   triangulation;
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> inv_mass_matrix;

    Vector<double>       density;
    Vector<double>       momentum;
    Vector<double>       energy;
    Vector<double>       density_old;
    Vector<double>       momentum_old;
    Vector<double>       energy_old;
    Vector<double>       rhs_density;
    Vector<double>       rhs_momentum;
    Vector<double>       rhs_energy;
   
    std::vector< Vector<double> > face_flux;
};

// Constructor
template <int dim>
NSProblem<dim>::NSProblem (unsigned int degree) :
    fe (degree),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
   
   xmin    = 0.0;
   xmax    = 1.0;
   n_cells = 500;

   dx      = (xmax - xmin) / n_cells;
   n_rk_stages = 3;
}

// Make grid and allocate memory for solution variables
template <int dim>
void NSProblem<dim>::make_grid_and_dofs ()
{
    GridGenerator::subdivided_hyper_cube (triangulation, n_cells, xmin, xmax);

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
    density_old.reinit (dof_handler.n_dofs());
    rhs_density.reinit (dof_handler.n_dofs());
   
    momentum.reinit (dof_handler.n_dofs());
    momentum_old.reinit (dof_handler.n_dofs());
    rhs_momentum.reinit (dof_handler.n_dofs());
   
    energy.reinit (dof_handler.n_dofs());
    energy_old.reinit (dof_handler.n_dofs());
    rhs_energy.reinit (dof_handler.n_dofs());   
   
    // Array to store flux across cell faces
    face_flux.resize(triangulation.n_active_cells()+1, Vector<double>(3));
}

// Set initial conditions
template <int dim>
void NSProblem<dim>::initialize ()
{
   std::cout << "Projecting initial condition ...\n";
   
   QGauss<dim>  quadrature_formula(fe.degree+1);
   
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   |
                            update_quadrature_points | 
                            update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   Vector<double>       cell_rhs_density  (dofs_per_cell);
   Vector<double>       cell_rhs_momentum (dofs_per_cell);
   Vector<double>       cell_rhs_energy   (dofs_per_cell);
   
   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   InitialCondition<dim> initial_condition;
   Vector<double> initial_value(3);
   double initial_density;
   double initial_momentum;
   double initial_energy;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      
      cell_rhs_density  = 0.0;
      cell_rhs_momentum = 0.0;
      cell_rhs_energy   = 0.0;
      
      
      // Flux integral over cell
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         // Get primitive variable at quadrature point
         initial_condition.vector_value(fe_values.quadrature_point(q_point),
                                        initial_value);
         // Convert primitive to conserved
         initial_density = initial_value(0);
         initial_momentum= initial_value(0) * initial_value(1);
         initial_energy  = initial_value(2)/(gas_gamma-1.0) + 
                           0.5 * initial_value(0) * pow(initial_value(1),2);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            cell_rhs_density(i) += (fe_values.shape_value (i, q_point) *
                                    initial_density *
                                    fe_values.JxW (q_point));
            cell_rhs_momentum(i)+= (fe_values.shape_value (i, q_point) *
                                    initial_momentum *
                                    fe_values.JxW (q_point));
            cell_rhs_energy(i)  += (fe_values.shape_value (i, q_point) *
                                    initial_energy *
                                    fe_values.JxW (q_point));
         }
      }
      
      
      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices (local_dof_indices);
      unsigned int ig, jg;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         ig = local_dof_indices[i];
         
         density (ig) = 0.0;
         momentum(ig) = 0.0;
         energy  (ig) = 0.0;
         
         for (unsigned int j=0; j<dofs_per_cell; ++j)
         {
            jg = local_dof_indices[j];
            density(ig)  += inv_mass_matrix(ig,jg) * cell_rhs_density(j);
            momentum(ig) += inv_mass_matrix(ig,jg) * cell_rhs_momentum(j);
            energy(ig)   += inv_mass_matrix(ig,jg) * cell_rhs_energy(j);
         }
         
      }
   }
}

// Assemble mass matrix for each cell
// Invert it and store
template <int dim>
void NSProblem<dim>::assemble_mass_matrix ()
{
    std::cout << "Constructing mass matrix ...\n";
    std::cout << "  Quadrature using " << fe.degree+1 << " points\n";

    QGauss<dim>  quadrature_formula(fe.degree+1);

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
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += fe_values.shape_value (i, q_point) *
                                        fe_values.shape_value (j, q_point) *
                                        fe_values.JxW (q_point);
       
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
void ns_flux (const double& density,
              const double& momentum,
              const double& energy,
              Vector<double>& flux)
{   
   double velocity = momentum / density;
   double pressure = (gas_gamma - 1.0) * (energy - 0.5 * momentum * velocity);
   flux(0) = momentum;
   flux(1) = pressure + momentum * velocity;
   flux(2) = (energy + pressure) * velocity;
}

// Lax-Friedrichs flux
void LaxFlux (Vector<double>& left_state,
              Vector<double>& right_state,
              Vector<double>& flux)
{
   Vector<double> left_flux(3);
   Vector<double> right_flux(3);
   
   double left_velocity = left_state(1) / left_state(0);
   double left_pressure = (gas_gamma-1.0) * (left_state(2) - 0.5 * left_state(1) * left_velocity );
   double left_sonic    = sqrt( gas_gamma * left_pressure / left_state(0) );
   double left_eig = fabs(left_velocity) + left_sonic;
   left_flux(0) = left_state(1);
   left_flux(1) = left_pressure + left_state(1) * left_velocity;
   left_flux(2) = (left_state(2) + left_pressure) * left_velocity;

   double right_velocity = right_state(1) / right_state(0);
   double right_pressure = (gas_gamma-1.0) * (right_state(2) - 0.5 * right_state(1) * right_velocity );
   double right_sonic    = sqrt( gas_gamma * right_pressure / right_state(0) );
   double right_eig = fabs(right_velocity) + right_sonic;
   right_flux(0) = right_state(1);
   right_flux(1) = right_pressure + right_state(1) * right_velocity;
   right_flux(2) = (right_state(2) + right_pressure) * right_velocity;
   
   
   double lambda = std::max ( left_eig, right_eig );
   
   for(unsigned int i=0; i<3; ++i)
   {
      flux(i) = 0.5 * ( left_flux(i) + right_flux(i) ) -
                0.5 * lambda * ( right_state(i) - left_state(i) );
   }
   
}

// Compute flux across cell faces
template <int dim>
void NSProblem<dim>::compute_face_flux ()
{
   const unsigned int dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> dofs(dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end(),
                                                  l_cell, r_cell;
   
   unsigned int l_dof, r_dof;
   Vector<double> left_state(3), right_state(3);
   
   // Loop over faces
   unsigned int n_faces = triangulation.n_active_cells() + 1;
   for (unsigned int i=0; i<n_faces; ++i)
   {
      if(i==0)
      {
         l_cell = cell;
         r_cell = cell;
         cell->get_dof_indices(dofs);
         l_dof = dofs[0];
         r_dof = dofs[0];
         
      }
      else if(i==n_faces-1)
      {
         l_cell = cell;
         r_cell = cell;
         cell->get_dof_indices(dofs);
         l_dof = dofs[dofs_per_cell-1];
         r_dof = dofs[dofs_per_cell-1];
         
      }
      else
      {
         l_cell = cell;
         r_cell = cell->neighbor(1);
         
         l_cell->get_dof_indices(dofs);
         l_dof = dofs[dofs_per_cell-1];
         
         r_cell->get_dof_indices(dofs);
         r_dof = dofs[0];
         
         ++cell;
      }
      
      left_state(0) = density  (l_dof);
      left_state(1) = momentum (l_dof);
      left_state(2) = energy   (l_dof);
      
      right_state(0) = density  (r_dof);
      right_state(1) = momentum (r_dof);
      right_state(2) = energy   (r_dof);
      
      LaxFlux (left_state, right_state, face_flux[i]);
   }
}

// Assemble system rhs
template <int dim>
void NSProblem<dim>::assemble_rhs ()
{
    QGaussLobatto<dim>  quadrature_formula(fe.degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | 
                             update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  density_values  (n_q_points);
    std::vector<double>  momentum_values (n_q_points);
    std::vector<double>  energy_values   (n_q_points);

    Vector<double>       cell_rhs_density  (dofs_per_cell);
    Vector<double>       cell_rhs_momentum (dofs_per_cell);
    Vector<double>       cell_rhs_energy   (dofs_per_cell);
   
    Vector<double>       flux(3);

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
            ns_flux(density_values[q_point], momentum_values[q_point], 
                    energy_values[q_point], flux);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                cell_rhs_density(i) += (fe_values.shape_grad (i, q_point)[0] *
                                        flux(0) *
                                        fe_values.JxW (q_point));
                cell_rhs_momentum(i)+= (fe_values.shape_grad (i, q_point)[0] *
                                        flux(1) *
                                        fe_values.JxW (q_point));
                cell_rhs_energy(i)  += (fe_values.shape_grad (i, q_point)[0] *
                                        flux(2) *
                                        fe_values.JxW (q_point));
            }
        }

        // Compute flux at cell boundaries TBD
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
           // Left face flux
           cell_rhs_density(i) += fe_values.shape_value (i, 0) *
                                  face_flux[cell->face_index(0)](0);
           cell_rhs_momentum(i)+= fe_values.shape_value (i, 0) *
                                  face_flux[cell->face_index(0)](1);
           cell_rhs_energy(i)  += fe_values.shape_value (i, 0) *
                                  face_flux[cell->face_index(0)](2);
           
           // Right face flux
           cell_rhs_density(i) -= fe_values.shape_value (i, n_q_points-1) *
                                  face_flux[cell->face_index(1)](0);
           cell_rhs_momentum(i)-= fe_values.shape_value (i, n_q_points-1) *
                                  face_flux[cell->face_index(1)](1);
           cell_rhs_energy(i)  -= fe_values.shape_value (i, n_q_points-1) *
                                  face_flux[cell->face_index(1)](2);
        }

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
void NSProblem<dim>::update (const unsigned int rk_stage)
{
   // Update conserved variables
   for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   {
      density(i)  = a_rk[rk_stage] * density_old(i) +
                    b_rk[rk_stage] * (density(i) + dt * rhs_density(i));
      momentum(i) = a_rk[rk_stage] * momentum_old(i) +
                    b_rk[rk_stage] * (momentum(i) + dt * rhs_momentum(i));
      energy(i)   = a_rk[rk_stage] * energy_old(i) +
                    b_rk[rk_stage] * (energy(i) + dt * rhs_energy(i));
   }

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
    assemble_mass_matrix ();
    initialize ();
    output_results ();

    dt = 0.1 * dx;
    double time = 0.0;
    unsigned int iter = 0;

    while (time < 0.2)
    {
       density_old  = density;
       momentum_old = momentum;
       energy_old   = energy;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
         compute_face_flux ();
         assemble_rhs ();
         update (rk);
       }
      time += dt;
      ++iter;
      std::cout << "Iter = " << iter << " time = " << time << endl;
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

