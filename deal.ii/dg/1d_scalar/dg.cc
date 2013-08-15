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
#include <numerics/vector_tools.h>
#include <numerics/matrix_tools.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>

#include <numerics/data_out.h>
#include <numerics/fe_field_function.h>

#include <base/convergence_table.h>

#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

// Number of variables: mass, momentum and energy
double xmin, xmax, xmid;
double u_left, u_right;
double beta = 1.0;
double final_time;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {lxf, kfvs};
enum TestCase {sine, sine2};

//------------------------------------------------------------------------------
// Scheme parameters
//------------------------------------------------------------------------------
struct Parameter
{
   int degree;
   double cfl;
   double final_time;
   TestCase test_case;
   unsigned int n_cells;
   unsigned int nstep;
};

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a, const double& b, const double& c)
{
   double result;
   if( a*b > 0.0 && b*c > 0.0)
   {
      result  = std::min( std::fabs(a), std::min(std::fabs(b), std::fabs(c)));
      result *= ((a>0.0) ? 1.0 : -1.0);
   }
   else 
   {
      result = 0.0;
   }
   
   return result;
}

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>() {}
   
   virtual void value (const Point<dim>   &p,
                       double& values) const;
};

// Initial condition
template<int dim>
void InitialCondition<dim>::value (const Point<dim>   &p,
                                   double& values) const
{
   double x = p[0];
   
   // test case: sine
   values = -std::sin (M_PI * x);
}

//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class Solution : public Function<dim>
{
public:
   Solution () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

// Exact solution
template<int dim>
double Solution<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
{
   double x = p[0] - final_time;
   
   // test case: sine
   return - std::sin (M_PI * x);

}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class ScalarProblem
{
public:
   ScalarProblem (Parameter param,
                  bool debug);
   void run ();
   
private:
   void solve ();
   void make_grid_and_dofs (unsigned int step);
   void initialize ();
   void assemble_mass_matrix ();
   void assemble_rhs ();
   void compute_averages ();
   void compute_dt ();
   void apply_limiter ();
   void update (const unsigned int rk_stage);
   void output_results (const double& time) const;
   void process_solution (unsigned int step);
   
   TestCase             test_case;
   bool                 debug;
   unsigned int         n_cells;
   double               dt;
   double               dx;
   double               cfl;
   double               min_residue;
   unsigned int         n_rk_stages;
   FluxType             flux_type;
   unsigned int         nstep;
   
   
   Triangulation<dim>   triangulation;
   FE_DGQ<dim>          fe;
   DoFHandler<dim>      dof_handler;
   
   SparsityPattern      sparsity_pattern;
   SparseMatrix<double> inv_mass_matrix;
   
   Vector<double>       solution;
   Vector<double>       solution_old;
   Vector<double>       rhs;
   
   std::vector< Vector<double> > average;
   
   double residual;
   double residual0;
   
   ConvergenceTable  convergence_table;
   typename DoFHandler<dim>::active_cell_iterator firstc, lastc;
   
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
ScalarProblem<dim>::ScalarProblem (Parameter param,
                                   bool debug) :
    test_case (param.test_case),
    debug(debug),
    n_cells (param.n_cells),
    cfl (param.cfl),
    nstep (param.nstep),
    fe (param.degree),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
   
   final_time = param.final_time;
   
   n_rk_stages = 3;
   flux_type = kfvs;
   
   if(test_case == sine)
   {
      xmin    = -1.0;
      xmax    = +1.0;
      min_residue= 1.0e20;
   }
   else if(test_case == sine2)
   {
      xmin    = -1.0;
      xmax    = +3.0;
      min_residue= 1.0e20;
   }
   else
   {
      std::cout << "Unknown test case\n";
   }
   
   dx = (xmax - xmin) / n_cells;

}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::make_grid_and_dofs (unsigned int step)
{
    if(step == 0)
       GridGenerator::subdivided_hyper_cube (triangulation, n_cells, xmin, xmax);
   else 
   {
      triangulation.refine_global (1);
      n_cells = triangulation.n_active_cells ();
      dx = (xmax - xmin)/n_cells;
   }

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
    solution.reinit (dof_handler.n_dofs());
    solution_old.reinit (dof_handler.n_dofs());
    rhs.reinit (dof_handler.n_dofs());
   
    average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
   
   // Find first and last cell
   // We need these for periodic bc
   // WARNING: This could be wrong with adaptive refinement.
   typename DoFHandler<dim>::active_cell_iterator
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();
   firstc = dof_handler.begin_active();
   for (unsigned int c=0; cell!=endc; ++cell, ++c)
   {
      if(c == triangulation.n_active_cells()-1)
         lastc = cell;
   }
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::initialize ()
{
   std::cout << "Projecting initial condition ...\n";
   
   QGauss<dim>  quadrature_formula(2 * fe.degree + 1);
   
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   |
                            update_quadrature_points | 
                            update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   Vector<double>       cell_rhs  (dofs_per_cell);
   
   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   InitialCondition<dim> initial_condition;
   double initial_value;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      
      cell_rhs  = 0.0;
      
      // Flux integral over cell
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         // Get primitive variable at quadrature point
         initial_condition.value(fe_values.quadrature_point(q_point),
                                        initial_value);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                            initial_value *
                            fe_values.JxW (q_point));
         }
      }
      
      
      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices (local_dof_indices);
      unsigned int ig, jg;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         ig = local_dof_indices[i];
         
         solution (ig) = 0.0;
         
         for (unsigned int j=0; j<dofs_per_cell; ++j)
         {
            jg = local_dof_indices[j];
            solution(ig)  += inv_mass_matrix(ig,jg) * cell_rhs(j);
         }
         
      }
   }
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::assemble_mass_matrix ()
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

//------------------------------------------------------------------------------
// Flux of the PDE model
//------------------------------------------------------------------------------
double physical_flux (const double& u)
{
   return u;
}

//------------------------------------------------------------------------------
// Lax-Friedrichs flux
//------------------------------------------------------------------------------
void LaxFlux (const double& left_state,
              const double& right_state,
              double& flux)
{
   flux = 0.0;
   cout << "Lax flux not finished\n";
   abort();
}

//------------------------------------------------------------------------------
// Convective KFVS split fluxes: sign=+1 give positive flux and
// sign=-1 gives negative flux
//------------------------------------------------------------------------------
void kfvs_c_split_flux (const double& u,
                        const int sign,
                        double& flux)
{
   double s, A, B;
   
   s    = sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   
   // inviscid flux
   flux = u * (A + B);
}

//------------------------------------------------------------------------------
// KFVS flux for navier-stokes
//------------------------------------------------------------------------------
void KFVSFlux (const double& left_state,
               const double& right_state,
               double& flux)
{
   double flux_c_pos;
   double flux_c_neg;
   
   kfvs_c_split_flux (left_state,  +1, flux_c_pos);
   kfvs_c_split_flux (right_state, -1, flux_c_neg);
   
  flux = flux_c_pos + flux_c_neg;
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void numerical_flux (const FluxType& flux_type,
                     double& left_state,
                     double& right_state,
                     double& flux)
{
   switch (flux_type) 
   {
      case lxf:
         LaxFlux (left_state, right_state, flux);
         break;
         
      case kfvs:
         KFVSFlux (left_state, right_state, flux);
         break;
         
      default:
         cout << "Unknown flux_type !!!\n";
         abort ();
   }
}

//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::assemble_rhs ()
{
    QGaussLobatto<dim>  quadrature_formula(fe.degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | 
                             update_JxW_values);

   // for getting neighbour cell solutions to compute intercell flux
   QTrapez<dim> quadrature_dummy;
   FEValues<dim> fe_values_neighbor (fe, quadrature_dummy,
                            update_values   | update_gradients);
   
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  solution_values  (n_q_points);
    std::vector< Tensor<1,dim> >  solution_grad (n_q_points);
   
   // for getting neighbor cell solution using trapezoidal rule
   std::vector<double>  solution_values_n  (2);   

    Vector<double>       cell_rhs  (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
    residual = 0.0;
   
    for (unsigned int c=0; cell!=endc; ++cell, ++c)
    {
        fe_values.reinit (cell);
       
        cell_rhs  = 0.0;

        // Compute conserved variables at quadrature points
        fe_values.get_function_values (solution,  solution_values);       

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
           double flux = physical_flux (solution_values[q_point]);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                cell_rhs(i) += (fe_values.shape_grad (i, q_point)[0] *
                                flux *
                                fe_values.JxW (q_point));
            }
        }
       
       // Computation of flux at cell boundaries
       double lf_left_state, lf_right_state;

        // left face flux
        // right state is from current cell
       lf_right_state = solution_values [0];
       
       if(cell == firstc)
       {
          // get left cell dof indices
          fe_values_neighbor.reinit (lastc);
          
          fe_values_neighbor.get_function_values (solution,  solution_values_n);
          
          lf_left_state = solution_values_n [1];
       }
       else
       {
          // get left cell dof indices
          fe_values_neighbor.reinit (cell->neighbor(0));
          
          fe_values_neighbor.get_function_values (solution,  solution_values_n);
          
          lf_left_state = solution_values_n [1];
       }
       
       double left_flux;
       numerical_flux (flux_type, lf_left_state, lf_right_state, 
                       left_flux);
              
       // right face flux
       // left state is from current cell
       double rf_left_state, rf_right_state;
       rf_left_state = solution_values [n_q_points-1];
       
       if(cell == lastc)
       {
          //right_state = u_right;
          fe_values_neighbor.reinit (firstc);
          
          fe_values_neighbor.get_function_values (solution,  solution_values_n);
          
          rf_right_state = solution_values_n [0];
       }
       else
       {          
          // get right cell to right face
          fe_values_neighbor.reinit (cell->neighbor(1));
          
          fe_values_neighbor.get_function_values (solution,  solution_values_n);
          
          rf_right_state = solution_values_n [0];
       }
       
       double right_flux;
       numerical_flux (flux_type, rf_left_state, rf_right_state, 
                       right_flux);

        // Add flux at cell boundaries
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
           // Left face flux
           cell_rhs(i) += fe_values.shape_value (i, 0) *
                          left_flux;
      
           
           // Right face flux
           cell_rhs(i) -= fe_values.shape_value (i, n_q_points-1) *
                          right_flux;
           
        }

        // Multiply by inverse mass matrix and add to rhs
        cell->get_dof_indices (local_dof_indices);
        unsigned int ig, jg;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            ig = local_dof_indices[i];
           
            rhs (ig) = 0.0;
           
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
               jg = local_dof_indices[j];
               rhs(ig)  += inv_mass_matrix(ig,jg) * cell_rhs(j);
            }
           
            residual += std::pow (rhs (ig), 2);
        }
       
    }

}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::compute_averages ()
{
   QGauss<dim>  quadrature_formula(fe.degree+1);
      
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | 
                            update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
      
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      
      average[c] = 0.0;
      for(unsigned int point=0; point<n_q_points; ++point)
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            average[c](0) += solution(local_dof_indices[i]) * 
                                     fe_values.shape_value (i, point) *
                                     fe_values.JxW (point);
         }
      
      average[c]  /= dx;
   }
}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::apply_limiter ()
{
   Assert (fe.degree==1, ExcIndexRange(fe.degree, 1, 2));
   
   const double bb = 1.5;
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   // dont limit in first cell, skip it.
   ++cell;
   
   double db, df, dl;
   for (unsigned int c=1; c<n_cells-1; ++c, ++cell)
   {
      cell->get_dof_indices (local_dof_indices);
      
      db = (average[c](0) - average[c-1](0)) / dx;
      df = (average[c+1](0) - average[c](0)) / dx;
      dl = minmod ( average[c](1), bb * db, bb * df);
      
      solution(local_dof_indices[0]) = average[c](0) - 
         0.5 * dx * dl;
      solution(local_dof_indices[1]) = average[c](0) + 
         0.5 * dx * dl;
      
   }
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::compute_dt ()
{
   dt = cfl * dx;
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::update (const unsigned int rk_stage)
{
   // Update conserved variables
   for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   {
      solution(i)  = a_rk[rk_stage] * solution_old(i) +
                     b_rk[rk_stage] * (solution(i) + dt * rhs(i));
   }

}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::output_results (const double& time) const
{
   static unsigned int c = 0;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");

    if(fe.degree <= 1)
       data_out.build_patches (1);
    else
       data_out.build_patches (2*fe.degree);
   
    std::string filename = "sol_" + Utilities::int_to_string(c) + ".gpl";
    ++c;
    cout << filename << endl;
   
   std::ofstream output (filename);
   data_out.write_gnuplot (output);

   // compute shear stress and heat flux at cell center
   // use midpoint quadtrature to get cell center values
   QMidpoint<dim> quadrature;
   FEValues<dim> fe_values (fe, quadrature, update_quadrature_points |
                            update_values   | update_gradients);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   std::vector<double> solution_values (1);
   
   Solution<dim> exact;
   double exact_solution;

   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   std::ofstream fo;
   fo.open ("cell_center.dat");
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, solution_values);
            
      Point<dim> x = fe_values.quadrature_point (0);
      exact_solution = exact.value(x);
      
      fo << fe_values.quadrature_point (0)(0) << " " 
         << exact_solution << "  "
         << solution_values[0] << endl;

   }
   
   fo.close ();
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::solve ()
{
    std::cout << "Solving 1-D scalar problem ...\n";

    //make_grid_and_dofs();
    assemble_mass_matrix ();
    initialize ();
    output_results (0.0);
    compute_averages ();

    double time = 0.0;
    unsigned int iter = 0;

    while (time < final_time)
    {
       solution_old  = solution;
       
       compute_dt ();
       if(time+dt > final_time) dt = final_time - time;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
         assemble_rhs ();
         update (rk);
         compute_averages ();
         //apply_limiter ();
       }
       
       if(iter==0)
       {
          std::cout << "Initial residual = " << residual << endl;
          residual0 = residual;
       }
       
       residual /= residual0;
       
      time += dt;
      ++iter;
       if(debug && iter % 10 == 0) output_results (time);
       
       if(debug)
      std::cout << "Iter = " << iter << " time = " << time 
                << " Res =" << residual << endl;
    }
   std::cout << "Iter = " << iter << " time = " << time 
   << " Res =" << residual << endl;

    output_results (time);
}
//------------------------------------------------------------------------------
// Compute error norms
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::process_solution (unsigned int step)
{
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(4),
                                      VectorTools::L2_norm);
   const double L2_error = difference_per_cell.l2_norm();
   
   const unsigned int n_active_cells=triangulation.n_active_cells();
   const unsigned int n_dofs=dof_handler.n_dofs();
   convergence_table.add_value("step", step);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("L2", L2_error);
   
   convergence_table.set_scientific("L2", true);
   
   std::cout << std::endl;
   convergence_table.write_text(std::cout);
}

//------------------------------------------------------------------------------
// Run with several refinements
// Compute error and convergence rate
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::run ()
{
   for(unsigned int step=0; step<nstep; ++step)
   {
      make_grid_and_dofs (step);
      solve ();
      process_solution (step);
   }

   convergence_table.set_precision("L2", 3);
   
   convergence_table.set_scientific("L2", true);
   
   convergence_table
   .evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
   
   convergence_table.set_tex_caption("cells", "\\# cells");
   convergence_table.set_tex_caption("dofs", "\\# dofs");
   convergence_table.set_tex_caption("L2", "L^2-error");
   
   convergence_table.set_tex_format("cells", "r");
   convergence_table.set_tex_format("dofs", "r");
   
   std::cout << std::endl;
   convergence_table.write_text(std::cout);
   
   std::ofstream error_table_file("error.tex");   
   convergence_table.write_tex(error_table_file);
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main ()
{
    deallog.depth_console (0);
    {
       Parameter param;
       param.degree = 1;
       param.n_cells = 100;
       param.nstep = 1;
       param.test_case = sine;
       param.cfl = 1.0/(2.0*param.degree+1.0);
       param.final_time = 10;
       
       bool debug = true;
       
       ScalarProblem<1> scalar_problem(param, debug);
       scalar_problem.run ();
    }

    return 0;
}

