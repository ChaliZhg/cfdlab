/* 1d DG code for euler equations. This is not particularly efficient code.
   TODO : Use MeshWorker to assemble rhs.
   * Legendre basis functions
   * TVD/TVB limiter
   * Momentum limiter (BDF)
   * Modified Moment limiter (BSB)
   * Characteristic based limiter
   * Positivity preserving limiter
   * KXRCF shock indicator
   * Numerical fluxes: Lax-Friedrich, KFVS
   *
   * Author: Praveen. C, http://praveen.cfdlab.net
*/
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_dgp.h>
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
#include <base/parameter_handler.h>
#include <base/convergence_table.h>

#include <numerics/data_out.h>
#include <numerics/fe_field_function.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

// Number of variables: mass, momentum and energy
const unsigned int n_var = 3;
double gas_gamma;
double gas_const;
double d_left, u_left, p_left;
double d_right, u_right, p_right;
double xmin, xmax, xmid;
double Mdx2; // for TVB limiter

// used in sedov problem
double pc, xc_l, xc_r;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
std::vector<double> a_rk, b_rk;

// Numerical flux functions
enum FluxType {lxf, kfvs, roe};
enum TestCase {sod, blast, blastwc, lax, shuosher, lowd, smooth};
enum ShockIndicator { ind_None, ind_density, ind_energy, ind_entropy };
enum ViscModel {constant, persson};

//------------------------------------------------------------------------------
// Computes entropy s = p / rho^gamma
//------------------------------------------------------------------------------
double entropy(double density, double momentum, double energy)
{
   double p = (gas_gamma-1.0) * (energy - 0.5 * momentum * momentum/density);
   return p / std::pow(density, gas_gamma);
}

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a, const double& b, const double& c)
{
   if(std::fabs(a) < Mdx2) return a;
   
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
// maxmod of two numbers
// Author: Sudarshan Kumar K
//------------------------------------------------------------------------------
double maxmod (const double& a, const double& b)
{   
   double result;
   if( a*b >= 0.0)
      result  = (std::fabs(a) > std::fabs(b)) ? a : b;
   else 
      result = 0.0;
   
   return result;
}

//------------------------------------------------------------------------------
// Compute matrix of eigenvectors = R
// and its inverse matrix = Ri
//------------------------------------------------------------------------------
void EigMat(double d, double m, double E, double R[][n_var], double Ri[][n_var])
{
   double v, p, c, h, f, g1, g2;

   g1 = gas_gamma - 1.0;
   g2 = g1 / 2.0;

   v = m / d;
   p = (gas_gamma - 1.0) * (E - 0.5 * d * v * v);
   c = sqrt(gas_gamma * p / d);
   h = c * c / g1 + 0.5 * v * v;
   f = d / c / 2.0;

   /* Inverse eigenvector-matrix */
   Ri[0][0] = 1.0 - g2 * v * v / c / c;
   Ri[1][0] = (g2 * v * v - v * c) / d / c;
   Ri[2][0] = -(g2 * v * v + v * c) / d / c;

   Ri[0][1] = g1 * v / c / c;
   Ri[1][1] = (c - g1 * v) / d / c;
   Ri[2][1] = (c + g1 * v) / d / c;

   Ri[0][2] = -g1 / c / c;
   Ri[1][2] = g1 / d / c;
   Ri[2][2] = -g1 / d / c;

   /* Eigenvector matrix */
   R[0][0] = 1.0;
   R[1][0] = v;
   R[2][0] = v * v / 2.0;

   R[0][1] = f;
   R[1][1] = (v + c) * f;
   R[2][1] = (h + v * c) * f;

   R[0][2] = -f;
   R[1][2] = -(v - c) * f;
   R[2][2] = -(h - v * c) * f;
}

//------------------------------------------------------------------------------
// U = R*U
//------------------------------------------------------------------------------
void Multi(double R[][n_var], std::vector<double>& U)
{
   std::vector<double> Ut(U);

   for(unsigned int i = 0; i < n_var; i++) 
   {
      U[i] = 0.0;
      for(unsigned int j = 0; j < n_var; j++)
         U[i] += R[i][j] * Ut[j];
   }
}
//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class ExactSolution : public Function<dim>
{
public:
   ExactSolution () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p, const unsigned int component = 0) const;
   virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                   const unsigned int  component = 0) const;
};

template<int dim>
double ExactSolution<dim>::value (const Point<dim>   &p, const unsigned int) const
{
   double x = p[0] - 2.0;
   return 1.0 + 0.1*sin(4.0*M_PI*x);
}

template <int dim>
Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim>   &p,
                                       const unsigned int) const
{
   Tensor<1,dim> return_value(0.0);
   double x = p[0] - 2.0;
   return_value[0] = 0.1 * (4.0*M_PI) * cos(4.0*M_PI*x);
   return return_value;
}
//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>() {}
   
   virtual void vector_value (const Point<dim>   &p,
                              Vector<double>& values) const;
   std::string test_case;
};

// Initial condition for density, velocity, pressure
template<int dim>
void InitialCondition<dim>::vector_value (const Point<dim>   &p,
                                          Vector<double>& values) const
{
   if(test_case == "sod" || 
      test_case == "lax" || 
      test_case == "lowd" || 
      test_case == "blastwc")
   {
      if(p[0] < xmid)
      {
         values(0) = d_left;
         values(1) = u_left;
         values(2) = p_left;
      }
      else
      {
         values(0) = d_right;
         values(1) = u_right;
         values(2) = p_right;
      }
   }
   else if(test_case == "sedov")
   {
      values(0) = 1.0;
      values(1) = 0.0;
      values(2) = 1.0e-12 * (gas_gamma-1.0);
      if(p[0] < xc_r && p[0] > xc_l)
         values(2) = pc;
   }
   else if(test_case == "blast")
   {
      if(p[0] < 0.1)
      {
         values(0) = 1;
         values(1) = 0;
         values(2) = 1000;
      }
      else if(p[0] > 0.9)
      {
         values(0) = 1;
         values(1) = 0;
         values(2) = 100;
      }
      else
      {
         values(0) = 1;
         values(1) = 0;
         values(2) = 0.01;
      }
   }
   else if(test_case == "shuosher")
   {
      if(p[0] < -4.0)
      {
         values(0) = 3.857143;
         values(1) = 2.629369;
         values(2) = 10.333333;
      }
      else
      {
         values(0) = 1 + 0.2*std::sin(5.0*p[0]);
         values(1) = 0;
         values(2) = 1;
      }
   }
   else if(test_case == "smooth")
   {
      values(0) = 1.0 + 0.1*sin(4.0*M_PI*p[0]);
      values(1) = 1.0;
      values(2) = 1.0;
   }
   else
   {
      std::cout << "Unknown test case\n";
      exit(0);
   }

}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class EulerProblem
{
public:
   EulerProblem (unsigned int degree, const ParameterHandler& prm);
   void run (double& h, int& ndof, double& L2_error, double& H1_error, double& Linf_error);
   
private:
   void make_grid_and_dofs ();
   void initialize ();
   void assemble_mass_matrix ();
   void compute_viscosity_constant ();
   void compute_viscosity_persson ();
   void compute_viscosity ();
   void assemble_rhs ();
   void compute_averages ();
   void compute_dt ();
   void identify_troubled_cells ();
   void apply_limiter ();
   void apply_limiter_TVB ();
   void apply_limiter_BDF ();
   void apply_limiter_BSB ();
   void apply_positivity_limiter ();
   void update (const unsigned int rk_stage);
   void output_results () const;
   void compute_errors (double& L2_error, double& H1_error, double& Linf_error) const;
   
   unsigned int         n_cells;
   std::string          test_case;
   double               dt;
   double               dx;
   double               cfl;
   double               cip;
   double               final_time;
   double               min_residue;
   unsigned int         max_iter;
   unsigned int         n_rk_stages;
   FluxType             flux_type;
   std::string          limiter;
   bool                 lim_char, lim_pos;
   ShockIndicator       shock_indicator;
   bool                 lbc_reflect, rbc_reflect, periodic;
   unsigned int         save_freq;
   ViscModel            visc_model;
   
   
   Triangulation<dim>   triangulation;
   FE_DGP<dim>          fe;
   DoFHandler<dim>      dof_handler;
   
   std::vector< Vector<double> > inv_mass_matrix;
   
   Vector<double>       density;
   Vector<double>       momentum;
   Vector<double>       energy;
   Vector<double>       density_old;
   Vector<double>       momentum_old;
   Vector<double>       energy_old;
   Vector<double>       rhs_density;
   Vector<double>       rhs_momentum;
   Vector<double>       rhs_energy;
   Vector<double>       viscosity;
   
   std::vector<double>  density_average;
   std::vector<double>  momentum_average;
   std::vector<double>  energy_average;
   
   std::vector<double>  residual;
   std::vector<double>  residual0;
   std::vector<bool>    is_troubled;
   unsigned int         n_troubled_cells;

   
   typename DoFHandler<dim>::active_cell_iterator firstc, lastc;
   std::vector<typename DoFHandler<dim>::active_cell_iterator> lcell, rcell;
   
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
EulerProblem<dim>::EulerProblem (unsigned int degree,
                                 const ParameterHandler& prm
                                 ) :
   fe (degree),
   dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));

   n_cells  = prm.get_integer("ncells");
   test_case= prm.get("test case");
   lim_char = prm.get_bool("characteristic limiter");
   lim_pos  = prm.get_bool("positivity limiter");
   cfl      = prm.get_double("cfl");
   cip      = prm.get_double("cip");
   double M = prm.get_double("M");
   save_freq= prm.get_integer("save frequency");
   limiter  = prm.get("limiter");
   std::string flux  = prm.get("flux");
   std::string indicator  = prm.get("indicator");
   max_iter = prm.get_integer("max iter");
   
   std::string vm  = prm.get("viscosity");
   if(vm == "constant") visc_model = constant;
   if(vm == "persson") visc_model = persson;
   
   if(limiter == "BDF" || limiter == "BSB") M = 0.0;
   
   n_rk_stages = std::min(degree,2u) + 1;
   a_rk.resize(n_rk_stages);
   b_rk.resize(n_rk_stages);
   if(n_rk_stages==1)
   {
      a_rk[0] = 0.0;
      b_rk[0] = 1.0;
   }
   else if(n_rk_stages==2)
   {
      a_rk = {0.0, 0.5};
      b_rk = {1.0, 0.5};
   }
   else if(n_rk_stages==3)
   {
      a_rk = {0.0, 3.0/4.0, 1.0/3.0};
      b_rk = {1.0, 1.0/4.0, 2.0/3.0};
   }
   else
   {
      std::cout << "This should not happen.\n";
      exit(0);
   }

   // Set flux enum type
   if(flux == "kfvs")
      flux_type = kfvs;
   else if(flux == "lxf")
      flux_type = lxf;
   else if(flux == "roe")
      flux_type = roe;

   if(indicator == "None")
      shock_indicator = ind_None;
   else if(indicator == "density")
      shock_indicator = ind_density;
   else if(indicator == "energy")
      shock_indicator = ind_energy;
   else if(indicator == "entropy")
      shock_indicator = ind_entropy;

   // If we are using shock indicator, then we should used TVD limiter
   if(shock_indicator != ind_None) M = 0.0;
   
   lbc_reflect = rbc_reflect = periodic = false;
   min_residue= 1.0e20;
   
   if(test_case == "sod")
   {
      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.5;
      final_time = 0.2;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      d_left  = 1.0;
      d_right = 0.125;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1.0;
      p_right = 0.1;
   }
   else if(test_case == "blast")
   {
      xmin    = 0.0;
      xmax    = 1.0;
      final_time = 0.038;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      lbc_reflect = rbc_reflect = true;
   }
   else if(test_case == "shuosher")
   {
      xmin    = -5;
      xmax    = 5;
      final_time = 1.8;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
   }
   else if(test_case == "lax")
   {
      xmin    = -5;
      xmax    = 5;
      xmid    = 0.0;
      final_time = 1.3;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      d_left  = 0.445;
      d_right = 0.5;
      
      u_left  = 0.698;
      u_right = 0.0;
      
      p_left  = 3.528;
      p_right = 0.571;
   }
   else if(test_case == "lowd")
   {
      xmin    = 0;
      xmax    = 1;
      xmid    = 0.5;
      final_time = 0.15;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = -2.0;
      u_right =  2.0;
      
      p_left  = 0.4;
      p_right = 0.4;
   }
   else if(test_case == "blastwc")
   {
      xmin    = 0;
      xmax    = 1.4;
      xmid    = 0.7;
      final_time = 0.012;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1000.0;
      p_right = 0.01;
   }
   else if(test_case == "sedov")
   {
      xmin    = -2.0;
      xmax    =  2.0;
      final_time = 0.001;
      
      gas_gamma = 1.4;
      gas_const = 1.0;

      AssertThrow( n_cells%2 == 1, ExcMessage("n_cells must be odd for sedov") );
      dx   = (xmax - xmin) / n_cells;
      // Energy in central cell
      double E0 = 3200000.0/dx;
      // Put in global variables
      pc   = E0 * (gas_gamma-1.0);
      xc_l = -0.5*dx;
      xc_r = +0.5*dx;
   }
   else if(test_case == "smooth")
   {
      xmin    = 0.0;
      xmax    = 1.0;
      final_time = 2.0;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      periodic = true;
   }
   else
   {
      std::cout << "Unknown test case\n";
   }
   
   cfl *= 1.0/(2.0*fe.degree+1.0);
   dx   = (xmax - xmin) / n_cells;
   Mdx2 = M * dx * dx;

}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::make_grid_and_dofs ()
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

   // allocate memory for inverse mass matrix. We store only diagonals.
   inv_mass_matrix.resize(n_cells);
   for (unsigned int c=0; c<n_cells; ++c)
      inv_mass_matrix[c].reinit(fe.degree+1);
   
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
   
    density_average.resize (triangulation.n_cells());
    momentum_average.resize (triangulation.n_cells());
    energy_average.resize (triangulation.n_cells());
    viscosity.reinit (triangulation.n_cells());
   
    residual.resize(3, 1.0);
    residual0.resize(3);
   
   is_troubled.resize(triangulation.n_cells());
   
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
   
   // for each cell find left cell and right cell
   lcell.resize(n_cells);
   rcell.resize(n_cells);
   cell = dof_handler.begin_active();
   for (unsigned int c=0; cell!=endc; ++cell, ++c)
   {
      if(c==0)
      {
         rcell[n_cells-1] = cell;
         lcell[c+1] = cell;
      }
      else if(c==n_cells-1)
      {
         rcell[c-1] = cell;
         lcell[0] = cell;
      }
      else
      {
         rcell[c-1] = cell;
         lcell[c+1] = cell;
      }
      
   }
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::initialize ()
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
   initial_condition.test_case = test_case;
   
   Vector<double> initial_value(n_var);
   double initial_density;
   double initial_momentum;
   double initial_energy;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (unsigned int c=0; cell!=endc; ++cell, ++c)
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
      unsigned int ig;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         ig = local_dof_indices[i];
         
         density(ig)  = inv_mass_matrix[c](i) * cell_rhs_density(i);
         momentum(ig) = inv_mass_matrix[c](i) * cell_rhs_momentum(i);
         energy(ig)   = inv_mass_matrix[c](i) * cell_rhs_energy(i);
      }
   }
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
// For Legendre basis, mass matrix is diagonal
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::assemble_mass_matrix ()
{
   std::cout << "Constructing mass matrix ...\n";
   std::cout << "  Quadrature using " << fe.degree+1 << " points\n";
   
   QGauss<dim>  quadrature_formula(fe.degree+1);
   
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values | update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   Vector<double>   cell_matrix (dofs_per_cell);
   
   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c=0; cell!=endc; ++cell,++c)
   {
      fe_values.reinit (cell);
      cell_matrix = 0.0;
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_matrix(i) += fe_values.shape_value (i, q_point) *
                              fe_values.shape_value (i, q_point) *
                              fe_values.JxW (q_point);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         inv_mass_matrix[c](i) = 1.0/cell_matrix(i);
   }
}


//------------------------------------------------------------------------------
// Flux for Euler equation
//------------------------------------------------------------------------------
void euler_flux (const double& density,
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

//------------------------------------------------------------------------------
// Lax-Friedrichs flux
//------------------------------------------------------------------------------
void LaxFlux (const Vector<double>& left_state,
              const Vector<double>& right_state,
              Vector<double>& flux)
{
   // Left state 
   double left_velocity = left_state(1) / left_state(0);
   double left_pressure = (gas_gamma-1.0) * (left_state(2) - 
                              0.5 * left_state(1) * left_velocity );
   double left_sonic    = sqrt( gas_gamma * left_pressure / left_state(0) );
   double left_eig      = fabs(left_velocity) + left_sonic;

   // Left flux
   Vector<double> left_flux(n_var);
   left_flux(0) = left_state(1);
   left_flux(1) = left_pressure + left_state(1) * left_velocity;
   left_flux(2) = (left_state(2) + left_pressure) * left_velocity;

   // Right state
   double right_velocity = right_state(1) / right_state(0);
   double right_pressure = (gas_gamma-1.0) * (right_state(2) - 
                              0.5 * right_state(1) * right_velocity );
   double right_sonic    = sqrt( gas_gamma * right_pressure / right_state(0) );
   double right_eig      = fabs(right_velocity) + right_sonic;

   // Right flux
   Vector<double> right_flux(n_var);
   right_flux(0) = right_state(1);
   right_flux(1) = right_pressure + right_state(1) * right_velocity;
   right_flux(2) = (right_state(2) + right_pressure) * right_velocity;
   
   // Maximum local wave speed at face
   double lambda = std::max ( left_eig, right_eig );
   
   for(unsigned int i=0; i<n_var; ++i)
      flux(i) = 0.5 * ( left_flux(i) + right_flux(i) ) -
                0.5 * lambda * ( right_state(i) - left_state(i) );
}

//------------------------------------------------------------------------------
// KFVS split fluxes: sign=+1 give positive flux and
// sign=-1 gives negative flux
// Copied from fv_ns_1d
//------------------------------------------------------------------------------
void kfvs_split_flux (const std::vector<double>& prim,
                      const int sign,
                      std::vector<double>& flux)
{
   double beta, s, A, B, E, fact;
   
   beta = 0.5 * prim[0] / prim[2];
   s    = prim[1] * sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   E    = prim[2]/(gas_gamma-1.0) + 0.5 * prim[0] * pow(prim[1], 2);
   fact = prim[1] * A + B;
   

   // inviscid flux
   flux[0] = prim[0] * fact;
   flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A +
             prim[0] * prim[1] * B;
   flux[2] = prim[1] * (E + prim[2]) * A +
             (E + 0.5 * prim[2]) * B;

}
//------------------------------------------------------------------------------
// KFVS flux
//------------------------------------------------------------------------------
void KFVSFlux (const Vector<double>& left_state,
               const Vector<double>& right_state,
               Vector<double>& flux)
{
   std::vector<double> left (n_var);
   std::vector<double> right (n_var);

   // Left primitive state 
   left[0] = left_state(0);
   left[1] = left_state(1) / left_state(0);
   left[2] = (gas_gamma-1.0) * (left_state(2) - 
                                    0.5 * left_state(1) * left[1] );
   
   // Right primitive state
   right[0] = right_state(0);
   right[1] = right_state(1) / right_state(0);
   right[2] = (gas_gamma-1.0) * (right_state(2) - 
                                    0.5 * right_state(1) * right[1] );
   
   std::vector<double> flux_pos (n_var);
   std::vector<double> flux_neg (n_var);
   
   kfvs_split_flux (left,  +1, flux_pos);
   kfvs_split_flux (right, -1, flux_neg);
   
   for(unsigned int i=0; i<n_var; ++i)
      flux(i) = flux_pos[i] + flux_neg[i];
}
//------------------------------------------------------------------------------
// Roe flux
//------------------------------------------------------------------------------
void ROEFlux (const Vector<double>& left_state,
              const Vector<double>& right_state,
              Vector<double>& flux)
{
   double left[n_var], right[n_var];
   
   left[0] = left_state[0];
   left[1] = left_state[1] / left_state[0];
   left[2] = (gas_gamma-1) * ( left_state[2] - 0.5 * left[0] * left[1] * left[1]);

   right[0] = right_state[0];
   right[1] = right_state[1] / right_state[0];
   right[2] = (gas_gamma-1) * ( right_state[2] - 0.5 * right[0] * right[1] * right[1]);
   
   double fl  = std::sqrt(left[0]);
   double fr  = std::sqrt(right[0]);
   double u   = (fl*left[1] + fr*right[1])/(fl + fr);
   
   double Hl  = gas_gamma*left[2]/left[0]/(gas_gamma-1.0) + 0.5*std::pow(left[1],2);
   double Hr  = gas_gamma*right[2]/right[0]/(gas_gamma-1.0) + 0.5*std::pow(right[1],2);
   
   // average of fluxes
   flux[0] = 0.5*(left[0]*left[1] + right[0]*right[1]);
   flux[1] = 0.5*(left[2] + left[0] * std::pow(left[1],2) +
                  right[2] + right[0] * std::pow(right[1],2));
   flux[2] = 0.5*(Hl*left[0]*left[1] + Hr*right[0]*right[1]);
   
   
   // Add conservative dissipation
   double H = (fl*Hl + fr*Hr)/(fl + fr);
   double a = std::sqrt((gas_gamma-1.0)*(H - 0.5*u*u));
   double R[3][3];
   R[0][0] = R[0][1] = R[0][2] = 1.0;
   R[1][0] = u-a; R[1][1] = u; R[1][2] = u + a;
   R[2][0] = H - u * a; R[2][1] = 0.5*u*u; R[2][2] = H + u * a;
   
   double Lambda[] = { fabs(u-a), fabs(u), fabs(u+a)};
   
   double dU[] = {
      right_state[0] - left_state[0],
      right_state[1] - left_state[1],
      right_state[2] - left_state[2]
   };
   
   double aa[3];
   aa[1] = (gas_gamma-1.0)/(a*a) * (dU[0]*(H-u*u) + u*dU[1] - dU[2]);
   aa[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a * aa[1]);
   aa[2] = dU[0] - aa[0] - aa[1];
   
   for(unsigned int i=0; i<3; ++i)
      for(unsigned int j=0; j<3; ++j)
         flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
}
//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void numerical_flux (const FluxType& flux_type,
                     Vector<double>& left_state,
                     Vector<double>& right_state,
                     Vector<double>& flux)
{
   switch (flux_type) 
   {
      case lxf: // Lax-Friedrich flux
         LaxFlux (left_state, right_state, flux);
         break;
         
      case kfvs: // Kinetic flux of Mandal & Deshpande
         KFVSFlux (left_state, right_state, flux);
         break;
         
      case roe: // Roe flux
         ROEFlux (left_state, right_state, flux);
         break;
         
      default:
         std::cout << "Unknown flux_type !!!\n";
         abort ();
   }
}


//------------------------------------------------------------------------------
//Convert conserved to primitive
//------------------------------------------------------------------------------
std::vector<double> con2prim (const Vector<double>& state)
{
   std::vector<double> prim (n_var);
   
   prim[0] = state(0);
   prim[1] = state(1) / state(0);
   prim[2] = (gas_gamma-1.0) * (state(2) - 
                                0.5 * state(1) * prim[1] );
   return prim;
}
//------------------------------------------------------------------------------
// constant artificial viscosity
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_viscosity_constant ()
{
   // constant viscosity per cell
   for(unsigned int i=0; i<n_cells; ++i)
   {
      double velocity = momentum_average[i] / density_average[i];
      double pressure = (gas_gamma-1.0) * ( energy_average[i] -
                                           0.5 * momentum_average[i] * velocity );
      double sonic = std::sqrt ( gas_gamma * pressure / density_average[i] );
      double speed = std::fabs(velocity) + sonic;
      viscosity(i) = dx * speed / fe.degree;
   }
}
//------------------------------------------------------------------------------
// Based on Persson and Peraire
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_viscosity_persson ()
{
   static const double s0 =  std::log10( 1.0/std::pow(fe.degree, 4) );
   static const double kappa = 4.0;
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   // constant viscosity per cell
   for(unsigned int c=0; cell!=endc; ++cell, ++c)
   {
      double velocity = momentum_average[c] / density_average[c];
      double pressure = (gas_gamma-1.0) * ( energy_average[c] -
                                           0.5 * momentum_average[c] * velocity );
      double sonic = std::sqrt ( gas_gamma * pressure / density_average[c] );
      double speed = std::fabs(velocity) + sonic;
      double mu0 = dx * speed / fe.degree;
      
      cell->get_dof_indices(dof_indices);
      double num = std::pow(density(dof_indices[fe.degree]),2);
      double den = 0;
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         den += std::pow(density(dof_indices[i]),2);
      
      double Se = num/den;
      double se = std::log10(Se);
      if(se < s0-kappa)
         viscosity(c) = 0;
      else if(se > s0+kappa)
         viscosity(c) = mu0;
      else
      {
         double arg = 0.5 * M_PI * (se - s0) / kappa;
         viscosity(c) = 0.5 * mu0 * (1.0 + std::sin(arg));
      }
   }
}
//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_viscosity ()
{
   if(limiter != "visc" || fe.degree == 0)
   {
      viscosity = 0;
      return;
   }
   
   switch(visc_model)
   {
      case constant:
         compute_viscosity_constant();
         break;
         
      case persson:
         compute_viscosity_persson();
         break;
   }
   
   /*
   Vector<double> tmp(n_cells);
   tmp = viscosity;
   for(unsigned int i=0; i<1; ++i)
   {
      for(unsigned int c=1; c<n_cells-1; ++c)
         viscosity(c) = 0.25*(tmp(c-1) + 2.0*tmp(c) + tmp(c+1));
   }
   */

}
//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::assemble_rhs ()
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

    std::vector<double>  density_values  (n_q_points);
    std::vector<double>  momentum_values (n_q_points);
    std::vector<double>  energy_values   (n_q_points);
   std::vector<Tensor<1,dim> >  density_grad  (n_q_points);
   std::vector<Tensor<1,dim> >  momentum_grad (n_q_points);
   std::vector<Tensor<1,dim> >  energy_grad   (n_q_points);
   
   // for getting neighbor cell solution using trapezoidal rule
   std::vector<double>  density_values_n  (2);
   std::vector<double>  momentum_values_n (2);
   std::vector<double>  energy_values_n   (2);
   
   std::vector<Tensor<1,dim> >  density_grad_n  (2);
   std::vector<Tensor<1,dim> >  momentum_grad_n (2);
   std::vector<Tensor<1,dim> >  energy_grad_n   (2);

    Vector<double>       cell_rhs_density  (dofs_per_cell);
    Vector<double>       cell_rhs_momentum (dofs_per_cell);
    Vector<double>       cell_rhs_energy   (dofs_per_cell);
   
    Vector<double>       flux(n_var);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
    residual[0] = residual[1] = residual[2] = 0.0;
   
    for (unsigned int c=0; cell!=endc; ++cell, ++c)
    {
        fe_values.reinit (cell);
       
        cell_rhs_density  = 0.0;
        cell_rhs_momentum = 0.0;
        cell_rhs_energy   = 0.0;

        // Compute conserved variables at quadrature points
        fe_values.get_function_values (density,  density_values);
        fe_values.get_function_values (momentum, momentum_values);
        fe_values.get_function_values (energy,   energy_values);
       
       fe_values.get_function_gradients (density,  density_grad);
       fe_values.get_function_gradients (momentum, momentum_grad);
       fe_values.get_function_gradients (energy,   energy_grad);

       
       // viscosity in this cell
       double mu = viscosity(c);

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            euler_flux(density_values[q_point], momentum_values[q_point],
                       energy_values[q_point], flux);
           flux(0) -= mu * density_grad[q_point][0];
           flux(1) -= mu * momentum_grad[q_point][0];
           flux(2) -= mu * energy_grad[q_point][0];
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
       
       // Computation of flux at cell boundaries
       Vector<double> lf_left_state(3), lf_right_state(3);
       Vector<double> lf_left_grad(3), lf_right_grad(3);
       double lf_mu_n;
       
        // left face flux
        // right state is from current cell
       lf_right_state(0) = density_values [0];
       lf_right_state(1) = momentum_values[0];
       lf_right_state(2) = energy_values  [0];
       
       lf_right_grad(0) = density_grad [0][0];
       lf_right_grad(1) = momentum_grad[0][0];
       lf_right_grad(2) = energy_grad  [0][0];
       
       if(c==0 && periodic==false)
       {
          if(lbc_reflect)
          {
             lf_left_state(0) = lf_right_state(0);
             lf_left_state(1) =-lf_right_state(1);
             lf_left_state(2) = lf_right_state(2);
          }
          else
          {
             lf_left_state(0) = lf_right_state(0);
             lf_left_state(1) = lf_right_state(1);
             lf_left_state(2) = lf_right_state(2);
          }
          lf_left_grad = lf_right_grad;
          lf_mu_n = mu;
       }
       else
       {
          // get left cell dof indices
          //fe_values_neighbor.reinit (cell->neighbor(0));
          fe_values_neighbor.reinit (lcell[c]);
          
          fe_values_neighbor.get_function_values (density,  density_values_n);
          fe_values_neighbor.get_function_values (momentum, momentum_values_n);
          fe_values_neighbor.get_function_values (energy,   energy_values_n);
          
          fe_values_neighbor.get_function_gradients (density,  density_grad_n);
          fe_values_neighbor.get_function_gradients (momentum, momentum_grad_n);
          fe_values_neighbor.get_function_gradients (energy,   energy_grad_n);
          
          lf_left_state(0) = density_values_n [1];
          lf_left_state(1) = momentum_values_n[1];
          lf_left_state(2) = energy_values_n  [1];
          
          lf_left_grad(0) = density_grad_n[1][0];
          lf_left_grad(1) = momentum_grad_n[1][0];
          lf_left_grad(2) = energy_grad_n[1][0];
          
          if(periodic)
             lf_mu_n = viscosity(n_cells-1);
          else
             lf_mu_n = viscosity(c-1);
       }
       
       Vector<double> left_flux(n_var);
       numerical_flux (flux_type, lf_left_state, lf_right_state, left_flux);
       double mu_l = 0.5*(mu + lf_mu_n);
       double cpen_l = cip * mu_l / dx;
       for(unsigned int i=0; i<n_var; ++i)
          left_flux(i) -= 0.5 * (lf_mu_n * lf_left_grad(i) + mu * lf_right_grad(i))
                          - cpen_l * (lf_right_state(i) - lf_left_state(i));
       
       // right face flux
       Vector<double> rf_left_state(n_var), rf_right_state(n_var);
       Vector<double> rf_left_grad(3), rf_right_grad(3);
       double rf_mu_n;

       // left state is from current cell
       rf_left_state(0) = density_values [n_q_points-1];
       rf_left_state(1) = momentum_values[n_q_points-1];
       rf_left_state(2) = energy_values  [n_q_points-1];
       
       rf_left_grad(0) = density_grad [n_q_points-1][0];
       rf_left_grad(1) = momentum_grad[n_q_points-1][0];
       rf_left_grad(2) = energy_grad  [n_q_points-1][0];
       
       if(c==triangulation.n_cells()-1 && periodic==false)
       {
          if(rbc_reflect)
          {
             rf_right_state(0) = rf_left_state(0);
             rf_right_state(1) =-rf_left_state(1);
             rf_right_state(2) = rf_left_state(2);
          }
          else
          {
             rf_right_state(0) = rf_left_state(0);
             rf_right_state(1) = rf_left_state(1);
             rf_right_state(2) = rf_left_state(2);
          }
          rf_right_grad = rf_left_grad;
          rf_mu_n = mu;
       }
       else
       {          
          // get right cell to right face
          //fe_values_neighbor.reinit (cell->neighbor(1));
          fe_values_neighbor.reinit (rcell[c]);
          
          fe_values_neighbor.get_function_values (density,  density_values_n);
          fe_values_neighbor.get_function_values (momentum, momentum_values_n);
          fe_values_neighbor.get_function_values (energy,   energy_values_n);
          
          fe_values_neighbor.get_function_gradients (density,  density_grad_n);
          fe_values_neighbor.get_function_gradients (momentum, momentum_grad_n);
          fe_values_neighbor.get_function_gradients (energy,   energy_grad_n);
          
          rf_right_state(0) = density_values_n [0];
          rf_right_state(1) = momentum_values_n[0];
          rf_right_state(2) = energy_values_n  [0];
          
          rf_right_grad(0) = density_grad_n [0][0];
          rf_right_grad(1) = momentum_grad_n[0][0];
          rf_right_grad(2) = energy_grad_n  [0][0];
          
          if(periodic)
             rf_mu_n = viscosity(0);
          else
             rf_mu_n = viscosity(c+1);
       }
       
       Vector<double> right_flux(3);
       numerical_flux (flux_type, rf_left_state, rf_right_state, right_flux);
       double mu_r = 0.5*(mu + rf_mu_n);
       double cpen_r = cip * mu_r / dx;
       for(unsigned int i=0; i<n_var; ++i)
          right_flux(i) -= 0.5 * (mu * rf_left_grad(i) + rf_mu_n * rf_right_grad(i))
                           - cpen_r * (rf_right_state(i) - rf_left_state(i));
       
        // Add flux at cell boundaries
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
           // Left face flux
           cell_rhs_density(i) += fe_values.shape_value (i, 0) *
                                  left_flux(0);
           cell_rhs_momentum(i)+= fe_values.shape_value (i, 0) *
                                  left_flux(1);
           cell_rhs_energy(i)  += fe_values.shape_value (i, 0) *
                                  left_flux(2);
           
           // Right face flux
           cell_rhs_density(i) -= fe_values.shape_value (i, n_q_points-1) *
                                  right_flux(0);
           cell_rhs_momentum(i)-= fe_values.shape_value (i, n_q_points-1) *
                                  right_flux(1);
           cell_rhs_energy(i)  -= fe_values.shape_value (i, n_q_points-1) *
                                  right_flux(2);
           
        }

        // Multiply by inverse mass matrix and add to rhs
        cell->get_dof_indices (local_dof_indices);
        unsigned int ig;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
           ig = local_dof_indices[i];
           
           rhs_density(ig)  = inv_mass_matrix[c](i) * cell_rhs_density(i);
           rhs_momentum(ig) = inv_mass_matrix[c](i) * cell_rhs_momentum(i);
           rhs_energy(ig)   = inv_mass_matrix[c](i) * cell_rhs_energy(i);
           
           residual[0] += std::pow (rhs_density (ig), 2);
           residual[1] += std::pow (rhs_momentum (ig), 2);
           residual[2] += std::pow (rhs_energy (ig), 2);
        }
       
    }

}

//------------------------------------------------------------------------------
// Compute cell average values
// For Legendre, first dof is cell average value
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_averages ()
{
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      cell->get_dof_indices (local_dof_indices);
      
      density_average[c]  = density  (local_dof_indices[0]);
      momentum_average[c] = momentum (local_dof_indices[0]);
      energy_average[c]   = energy   (local_dof_indices[0]);
   }
}

//------------------------------------------------------------------------------
// KXRCF shock indicator
// You can choose density or entropy as the indicator variables
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::identify_troubled_cells ()
{
   // If no shock indicator, then we apply limiter in every cell
   if(shock_indicator == ind_None)
   {
      for(unsigned int c=0; c<n_cells; ++c)
         is_troubled[c] = true;
      return;
   }

   QTrapez<dim>  quadrature_formula;

   FEValues<dim> fe_values (fe, quadrature_formula, update_values);
   FEValues<dim> fe_values_nbr (fe, quadrature_formula, update_values);
   std::vector<double> density_face_values(2), momentum_face_values(2),
                       energy_face_values(2);
   std::vector<double> density_face_values_n(2), momentum_face_values_n(2),
                       energy_face_values_n(2);
   
   const double dxfactor = std::pow( dx, 0.5*(fe.degree+1) );
   double density_nbr_l, momentum_nbr_l, energy_nbr_l;
   double density_nbr_r, momentum_nbr_r, energy_nbr_r;
   double ent_l, ent_r, ent_nbr_l, ent_nbr_r;
   double var_avg;

   n_troubled_cells = 0;

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(density, density_face_values);
      fe_values.get_function_values(momentum, momentum_face_values);
      fe_values.get_function_values(energy, energy_face_values);
      
      // Use cell average velocity to decide inflow/outflow boundary
      double u = momentum_average[c] / density_average[c];
      
      if(u > 0.0) // left face is inflow face
      {
         if(c==0 && periodic==false)
         {
            if(lbc_reflect)
            {
               density_nbr_l = density_face_values[0];
               momentum_nbr_l = -momentum_face_values[0];
               energy_nbr_l = energy_face_values[0];
            }
            else
            {
               density_nbr_l = density_face_values[0];
               momentum_nbr_l = momentum_face_values[0];
               energy_nbr_l = energy_face_values[0];
            }
         }
         else
         {
            fe_values_nbr.reinit (lcell[c]);
          
            fe_values_nbr.get_function_values (density,  density_face_values_n);
            fe_values_nbr.get_function_values (momentum, momentum_face_values_n);
            fe_values_nbr.get_function_values (energy,   energy_face_values_n);

            density_nbr_l = density_face_values_n[1];
            momentum_nbr_l = momentum_face_values_n[1];
            energy_nbr_l = energy_face_values_n[1];
         }
      }
      
      if(u < 0.0) // right face is inflow face
      {
         if(c==n_cells-1 && periodic==false)
         {
            if(rbc_reflect)
            {
               density_nbr_r = density_face_values[1];
               momentum_nbr_r = -momentum_face_values[1];
               energy_nbr_r = energy_face_values[1];
            }
            else
            {
               density_nbr_r = density_face_values[1];
               momentum_nbr_r = momentum_face_values[1];
               energy_nbr_r = energy_face_values[1];
            }
         }
         else
         {
            fe_values_nbr.reinit (rcell[c]);
          
            fe_values_nbr.get_function_values (density,  density_face_values_n);
            fe_values_nbr.get_function_values (momentum, momentum_face_values_n);
            fe_values_nbr.get_function_values (energy,   energy_face_values_n);

            density_nbr_r = density_face_values_n[0];
            momentum_nbr_r = momentum_face_values_n[0];
            energy_nbr_r = energy_face_values_n[0];
         }

      }

      double ind = 0;

      switch(shock_indicator)
      {
         case ind_density:
            ind = (u > 0.0) ? density_face_values[0] - density_nbr_l : 0.0
                + (u < 0.0) ? density_face_values[1] - density_nbr_r : 0.0;
            var_avg = density_average[c];
            break;
         case ind_energy:
            ind = (u > 0.0) ? energy_face_values[0] - energy_nbr_l : 0.0
                + (u < 0.0) ? energy_face_values[1] - energy_nbr_r : 0.0;
            var_avg = energy_average[c];
            break;
         case ind_entropy:
            ent_l = entropy(density_face_values[0], momentum_face_values[0], energy_face_values[0]);
            ent_r = entropy(density_face_values[1], momentum_face_values[1], energy_face_values[1]);
            ent_nbr_l = entropy(density_nbr_l, momentum_nbr_l, energy_nbr_l);
            ent_nbr_r = entropy(density_nbr_r, momentum_nbr_r, energy_nbr_r);
            ind = (u > 0.0) ? ent_l - ent_nbr_l : 0.0
                + (u < 0.0) ? ent_r - ent_nbr_r : 0.0;
            var_avg = entropy(density_average[c], momentum_average[c], energy_average[c]);
            break;
         default:
            exit(0);
      }

      ind = std::fabs(ind / (dxfactor * var_avg) );
      is_troubled[c] = false;
      if(ind > 1.0)
      {
         is_troubled[c] = true;
         ++n_troubled_cells;
      }
   }

   std::cout << "Number of troubled cells = " << n_troubled_cells << std::endl;
}

//------------------------------------------------------------------------------
// Apply chosen limiter
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_limiter ()
{
   if(limiter == "TVB")
      apply_limiter_TVB ();
   else if(limiter == "BDF")
      apply_limiter_BDF ();
   else if(limiter == "BSB")
      apply_limiter_BSB ();
   else if(limiter == "None" || limiter == "visc")
      return;
   else
   {
      std::cout << "Unknown limiter\n";
      exit(0);
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_limiter_TVB ()
{
   if(fe.degree == 0) return;
   
   QTrapez<dim>  quadrature_formula;
   
   FEValues<dim> fe_values (fe, quadrature_formula, update_values);
   std::vector<double> density_face_values(2), momentum_face_values(2),
                       energy_face_values(2);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   std::vector<double> db(n_var), df(n_var), DB(n_var), DF(n_var);
   std::vector<double> dl(n_var), dr(n_var);
   double density_left, density_right;
   double momentum_left, momentum_right;
   double energy_left, energy_right;
   
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   if(is_troubled[c])
   {
      fe_values.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      fe_values.get_function_values(density, density_face_values);
      fe_values.get_function_values(momentum, momentum_face_values);
      fe_values.get_function_values(energy, energy_face_values);
      
      unsigned int lc = (c==0) ? n_cells-1 : c-1;
      unsigned int rc = (c==n_cells-1) ? 0 : c+1;
      
      if(c==0 && !periodic)
      {
         density_left = density_average[c];
         if(lbc_reflect)
            momentum_left = -momentum_average[c];
         else
            momentum_left = momentum_average[c];
         energy_left = energy_average[c];
         
         density_right = density_average[c+1];
         momentum_right = momentum_average[c+1];
         energy_right = energy_average[c+1];
      }
      else if(c == n_cells-1 && !periodic)
      {
         density_left = density_average[c-1];
         momentum_left = momentum_average[c-1];
         energy_left = energy_average[c-1];
         
         density_right = density_average[c];
         if(rbc_reflect)
            momentum_right = -momentum_average[c];
         else
            momentum_right = momentum_average[c];
         energy_right = energy_average[c];
      }
      else
      {
         density_left = density_average[lc];
         momentum_left = momentum_average[lc];
         energy_left = energy_average[lc];
         
         density_right = density_average[rc];
         momentum_right = momentum_average[rc];
         energy_right = energy_average[rc];
      }
      
      // density
      db[0] = density_average[c] - density_left;
      df[0] = density_right - density_average[c];
      DB[0] = density_average[c] - density_face_values[0];
      DF[0] = density_face_values[1] - density_average[c];
      
      // momentum
      db[1] = momentum_average[c] - momentum_left;
      df[1] = momentum_right - momentum_average[c];
      DB[1] = momentum_average[c] - momentum_face_values[0];
      DF[1] = momentum_face_values[1] - momentum_average[c];
      
      // energy
      db[2] = energy_average[c] - energy_left;
      df[2] = energy_right - energy_average[c];
      DB[2] = energy_average[c] - energy_face_values[0];
      DF[2] = energy_face_values[1] - energy_average[c];

      double R[n_var][n_var], Ri[n_var][n_var];
      if(lim_char)
      {
         EigMat(density_average[c], 
                momentum_average[c], 
                energy_average[c], R, Ri);
         Multi(Ri, db);
         Multi(Ri, df);
         Multi(Ri, DB);
         Multi(Ri, DF);
      }

      double diff = 0;
      for(unsigned int i=0; i<n_var; ++i)
      {
         dl[i] = minmod (DB[i], db[i], df[i]);
         dr[i] = minmod (DF[i], db[i], df[i]);
         diff += std::fabs(dl[i] - DB[i]) + std::fabs(dr[i]-DF[i]);
      }
      diff /= (2*n_var);

      // If diff is nonzero, then limiter is active.
      // Then we keep only linear part
      if(diff > 1.0e-10)
      {
         if(lim_char) 
         {
            Multi(R, dl);
            Multi(R, dr);
         }
         density(local_dof_indices[1])  = 0.5*(dl[0] + dr[0]) / fe_values.shape_value(1,1);
         momentum(local_dof_indices[1]) = 0.5*(dl[1] + dr[1]) / fe_values.shape_value(1,1);
         energy(local_dof_indices[1])   = 0.5*(dl[2] + dr[2]) / fe_values.shape_value(1,1);
         // Higher dofs are set to zero
         for(unsigned int i=2; i<dofs_per_cell; ++i)
         {
            density(local_dof_indices[i])  = 0.0;
            momentum(local_dof_indices[i]) = 0.0;
            energy(local_dof_indices[i])   = 0.0;
         }
      }
      
   }
}
//------------------------------------------------------------------------------
// Apply moment limiter of Biswas, Devine, Flaherty
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_limiter_BDF ()
{
   if(fe.degree == 0) return;
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   std::vector<unsigned int> left_dof_indices (dofs_per_cell);
   std::vector<unsigned int> right_dof_indices (dofs_per_cell);
   
   std::vector< std::vector<double> > db(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > df(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > DC(dofs_per_cell, std::vector<double>(n_var));

   // Temporary storage
   Vector<double> density_n(density);
   Vector<double> momentum_n(momentum);
   Vector<double> energy_n(energy);

   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   if(is_troubled[c])
   {
      cell->get_dof_indices (local_dof_indices);
      
      if(c==0 && !periodic)
      {
         rcell[c]->get_dof_indices (right_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = 0.0;
            if(lbc_reflect)
               db[i][1] = 2.0*momentum(local_dof_indices[i]);
            else
               db[i][1] = 0.0;
            db[i][2] = 0.0;

            df[i][0] = density(right_dof_indices[i]) - density(local_dof_indices[i]);
            df[i][1] = momentum(right_dof_indices[i]) - momentum(local_dof_indices[i]);
            df[i][2] = energy(right_dof_indices[i]) - energy(local_dof_indices[i]);
         }
      }
      else if(c == n_cells-1 && !periodic)
      {
         lcell[c]->get_dof_indices (left_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = density(local_dof_indices[i]) - density(left_dof_indices[i]);
            db[i][1] = momentum(local_dof_indices[i]) - momentum(left_dof_indices[i]);
            db[i][2] = energy(local_dof_indices[i]) - energy(left_dof_indices[i]);

            df[i][0] = 0.0;
            if(rbc_reflect)
               df[i][1] = 2.0*momentum(local_dof_indices[i]);
            else
               df[i][1] = 0.0;
            df[i][2] = 0.0;
         }
      }
      else
      {
         lcell[c]->get_dof_indices (left_dof_indices);
         rcell[c]->get_dof_indices (right_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = density(local_dof_indices[i]) - density(left_dof_indices[i]);
            db[i][1] = momentum(local_dof_indices[i]) - momentum(left_dof_indices[i]);
            db[i][2] = energy(local_dof_indices[i]) - energy(left_dof_indices[i]);

            df[i][0] = density(right_dof_indices[i]) - density(local_dof_indices[i]);
            df[i][1] = momentum(right_dof_indices[i]) - momentum(local_dof_indices[i]);
            df[i][2] = energy(right_dof_indices[i]) - energy(local_dof_indices[i]);
         }
      }
      
      for(unsigned int i=0; i<dofs_per_cell; ++i)
      {
         DC[i][0] = density(local_dof_indices[i]);
         DC[i][1] = momentum(local_dof_indices[i]);
         DC[i][2] = energy(local_dof_indices[i]);
      }

      double R[n_var][n_var], Ri[n_var][n_var];
      if(lim_char)
      {
         EigMat(density_average[c], 
                momentum_average[c], 
                energy_average[c], R, Ri);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            Multi(Ri, db[i]);
            Multi(Ri, df[i]);
            Multi(Ri, DC[i]);
         }
      }

      // Legendre in deal.ii is normalized. Moment limiter is BDF paper is
      // given for non-normalized basis functions. We apply correct
      // transformation here to account for this difference
      // cell average value is unchanged
      bool to_limit = true;
      for(unsigned int i=dofs_per_cell-1; i>=1; --i)
      {
         if(to_limit)
         {
            double l = (2*i - 1)*std::sqrt(2*i+1);
            double s = std::sqrt(2*i-1);
            std::vector<double> dcn(n_var);
            dcn[0] = minmod(l*DC[i][0], s*db[i-1][0], s*df[i-1][0])/l;
            dcn[1] = minmod(l*DC[i][1], s*db[i-1][1], s*df[i-1][1])/l;
            dcn[2] = minmod(l*DC[i][2], s*db[i-1][2], s*df[i-1][2])/l;
            double diff = std::fabs(dcn[0]-DC[i][0]) 
                        + std::fabs(dcn[1]-DC[i][1])
                        + std::fabs(dcn[2]-DC[i][2]);
            if(lim_char) Multi(R, dcn);
            density_n(local_dof_indices[i])  = dcn[0];
            momentum_n(local_dof_indices[i]) = dcn[1];
            energy_n(local_dof_indices[i])   = dcn[2];
            if(diff < 1.0e-10) to_limit = false; // Remaining dofs will not be limited
         }
         else
         {
            density_n(local_dof_indices[i])  = density(local_dof_indices[i]); 
            momentum_n(local_dof_indices[i]) = momentum(local_dof_indices[i]);
            energy_n(local_dof_indices[i])   = energy(local_dof_indices[i]);
         }
      }

   }

   // Now copy to main arrays
   density = density_n;
   momentum= momentum_n;
   energy  = energy_n;
}

//------------------------------------------------------------------------------
// Moment limiter of Burbeau, Sagaut, Bruneau
// Author: Sudarshan Kumar K
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_limiter_BSB ()
{
   if(fe.degree == 0) return;
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   std::vector<unsigned int> left_dof_indices (dofs_per_cell);
   std::vector<unsigned int> right_dof_indices (dofs_per_cell);
   
   std::vector< std::vector<double> > db(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > df(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > DC(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > DC_l(dofs_per_cell, std::vector<double>(n_var));
   std::vector< std::vector<double> > DC_r(dofs_per_cell, std::vector<double>(n_var));

   // Temporary storage
   Vector<double> density_n(density);
   Vector<double> momentum_n(momentum);
   Vector<double> energy_n(energy);

   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   if(is_troubled[c])
   {
      cell->get_dof_indices (local_dof_indices);
      
      if(c==0 && !periodic)
      {
         rcell[c]->get_dof_indices (right_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = 0.0;
            DC_l[i][0]=density(local_dof_indices[i]);
            if(lbc_reflect)
            {
               db[i][1] = 2.0*momentum(local_dof_indices[i]);
               DC_l[i][1]=-momentum(local_dof_indices[i]);
            }
            else
            {
               db[i][1] = 0.0;
               DC_l[i][1]=momentum(local_dof_indices[i]);
            }
            db[i][2] = 0.0;
            DC_l[i][2]=energy(local_dof_indices[i]);
            
            df[i][0] = density(right_dof_indices[i]) - density(local_dof_indices[i]);
            df[i][1] = momentum(right_dof_indices[i]) - momentum(local_dof_indices[i]);
            df[i][2] = energy(right_dof_indices[i]) - energy(local_dof_indices[i]);
            
            DC_r[i][0] = density(right_dof_indices[i]);
            DC_r[i][1] = momentum(right_dof_indices[i]);
            DC_r[i][2] = energy(right_dof_indices[i]);
         }
      }
      else if(c == n_cells-1 && !periodic)
      {
         lcell[c]->get_dof_indices (left_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = density(local_dof_indices[i]) - density(left_dof_indices[i]);
            db[i][1] = momentum(local_dof_indices[i]) - momentum(left_dof_indices[i]);
            db[i][2] = energy(local_dof_indices[i]) - energy(left_dof_indices[i]);
            
            DC_l[i][0] = density(left_dof_indices[i]);
            DC_l[i][1] = momentum(left_dof_indices[i]);
            DC_l[i][2] = energy(left_dof_indices[i]);

            df[i][0] = 0.0;
            DC_r[i][0]=density(local_dof_indices[i]);
            if(rbc_reflect)
            {
               df[i][1] = 2.0*momentum(local_dof_indices[i]);
               DC_r[i][1]=-momentum(local_dof_indices[i]);
            }
            else
            {
               df[i][1] = 0.0;
               DC_r[i][1]=momentum(local_dof_indices[i]);
            }
            df[i][2] = 0.0;
            DC_r[i][2]=energy(local_dof_indices[i]);
         }
      }
      else
      {
         lcell[c]->get_dof_indices (left_dof_indices);
         rcell[c]->get_dof_indices (right_dof_indices);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            db[i][0] = density(local_dof_indices[i]) - density(left_dof_indices[i]);
            db[i][1] = momentum(local_dof_indices[i]) - momentum(left_dof_indices[i]);
            db[i][2] = energy(local_dof_indices[i]) - energy(left_dof_indices[i]);
            
            df[i][0] = density(right_dof_indices[i]) - density(local_dof_indices[i]);
            df[i][1] = momentum(right_dof_indices[i]) - momentum(local_dof_indices[i]);
            df[i][2] = energy(right_dof_indices[i]) - energy(local_dof_indices[i]);
            
            DC_l[i][0] = density(left_dof_indices[i]);
            DC_l[i][1] = momentum(left_dof_indices[i]);
            DC_l[i][2] = energy(left_dof_indices[i]);
            
            DC_r[i][0] = density(right_dof_indices[i]);
            DC_r[i][1] = momentum(right_dof_indices[i]);
            DC_r[i][2] = energy(right_dof_indices[i]);
         }
      }
      
      for(unsigned int i=0; i<dofs_per_cell; ++i)
      {
         DC[i][0] = density(local_dof_indices[i]);
         DC[i][1] = momentum(local_dof_indices[i]);
         DC[i][2] = energy(local_dof_indices[i]);
      }

      double R[n_var][n_var], Ri[n_var][n_var];
      if(lim_char)
      {
         EigMat(density_average[c], 
                momentum_average[c], 
                energy_average[c], R, Ri);
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            Multi(Ri, db[i]);
            Multi(Ri, df[i]);
            Multi(Ri, DC[i]);
            Multi(Ri, DC_l[i]);
            Multi(Ri, DC_r[i]);
         }
      }

      // Legendre in deal.ii is normalized. Moment limiter is BDF & BSB paper is
      // given for non-normalized basis functions. We apply correct
      // transformation here to account for this difference
      // cell average value is unchanged
      bool to_limit = true;
      for(unsigned int i=dofs_per_cell-1; i>=1; --i)
      {
         if(to_limit)
         {
            double l = (2*i - 1)*std::sqrt(2*i+1);
            double s = std::sqrt(2*i-1);
            double t = std::sqrt(2*i+1);
            std::vector<double> dcn(n_var),dcn_m(n_var);
            std::vector<double> ur(n_var),ul(n_var),u_max(n_var);
            dcn[0] = minmod(l*DC[i][0], s*db[i-1][0], s*df[i-1][0])/l;
            dcn[1] = minmod(l*DC[i][1], s*db[i-1][1], s*df[i-1][1])/l;
            dcn[2] = minmod(l*DC[i][2], s*db[i-1][2], s*df[i-1][2])/l;

           double diff = std::fabs(dcn[0]-DC[i][0]) 
                       + std::fabs(dcn[1]-DC[i][1])
                       + std::fabs(dcn[2]-DC[i][2]);
            if(diff>1.0e-10)
            {
               ur[0] = s*DC_r[i-1][0] - (2*i-1)*t*DC_r[i][0];
               ur[1] = s*DC_r[i-1][1] - (2*i-1)*t*DC_r[i][1];
               ur[2] = s*DC_r[i-1][2] - (2*i-1)*t*DC_r[i][2];
               
               ul[0] = s*DC_l[i-1][0] + (2*i-1)*t*DC_l[i][0];
               ul[1] = s*DC_l[i-1][1] + (2*i-1)*t*DC_l[i][1];
               ul[2] = s*DC_l[i-1][2] + (2*i-1)*t*DC_l[i][2];
               
               u_max[0] = minmod(l*DC[i][0], (ur[0]-s*DC[i-1][0]), (s*DC[i-1][0]-ul[0]) ) /l;
               u_max[1] = minmod(l*DC[i][1], (ur[1]-s*DC[i-1][1]), (s*DC[i-1][1]-ul[1]) ) /l;
               u_max[2] = minmod(l*DC[i][2], (ur[2]-s*DC[i-1][2]), (s*DC[i-1][2]-ul[2]) ) /l;
               
               dcn_m[0] = maxmod(dcn[0],u_max[0]);
               dcn_m[1] = maxmod(dcn[1],u_max[1]);
               dcn_m[2] = maxmod(dcn[2],u_max[2]);
               if(lim_char) Multi(R, dcn_m);
               density_n(local_dof_indices[i])  = dcn_m[0];
               momentum_n(local_dof_indices[i]) = dcn_m[1];
               energy_n(local_dof_indices[i])   = dcn_m[2];
            }
            else
            {
               density_n(local_dof_indices[i])  = density(local_dof_indices[i]); 
               momentum_n(local_dof_indices[i]) = momentum(local_dof_indices[i]);
               energy_n(local_dof_indices[i])   = energy(local_dof_indices[i]);
               to_limit = false; // Remaining dofs will not be limited
            }
         }
         else
         {
            density_n(local_dof_indices[i])  = density(local_dof_indices[i]); 
            momentum_n(local_dof_indices[i]) = momentum(local_dof_indices[i]);
            energy_n(local_dof_indices[i])   = energy(local_dof_indices[i]);
         }
      }

   }

   // Now copy to main arrays
   density = density_n;
   momentum= momentum_n;
   energy  = energy_n;
}

//------------------------------------------------------------------------------
// Apply positivity limiter
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_positivity_limiter ()
{
   if(fe.degree == 0) return;
   
   // Need 2N - 3 >= degree for the quadrature to be exact.
   unsigned int N = (fe.degree + 3)/2;
   if((fe.degree+3)%2 != 0) N += 1;
   QGaussLobatto<dim>  quadrature_formula(N);
   const unsigned int n_q_points = quadrature_formula.size();
   FEValues<dim> fe_values (fe, quadrature_formula, update_values);
   std::vector<double> density_values(n_q_points), momentum_values(n_q_points),
                       energy_values(n_q_points);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   double eps = 1.0e-13;
   for (unsigned int c=0; c<n_cells; ++c)
   {
      double velocity = momentum_average[c] / density_average[c];
      double pressure = (gas_gamma-1.0) * ( energy_average[c] -
                                           0.5 * momentum_average[c] * velocity );
      eps = std::min(eps, density_average[c]);
      eps = std::min(eps, pressure);
   }
   if(eps < 0.0)
   {
      std::cout << "Fatal: Negative states\n";
      exit(0);
   }

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   for (unsigned int c=0; cell!=endc; ++cell, ++c)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      
      // First limit density
      fe_values.get_function_values(density, density_values);
      
      // find minimum density at GLL points
      double rho_min = 1.0e20;
      for(unsigned int q=0; q<n_q_points; ++q)
         rho_min = std::min(rho_min, density_values[q]);
      
      double rat = std::fabs(density_average[c] - eps) /
                   (std::fabs(density_average[c] - rho_min) + 1.0e-13);
      double theta1 = std::min(rat, 1.0);
      
      for(unsigned int i=1; i<dofs_per_cell; ++i)
         density(local_dof_indices[i]) *= theta1;
      
      // now limit pressure
      fe_values.get_function_values(density, density_values);
      fe_values.get_function_values(momentum, momentum_values);
      fe_values.get_function_values(energy, energy_values);
      
      double theta2 = 1.0;
      for(unsigned int q=0; q<n_q_points; ++q)
      {
         double pressure = (gas_gamma-1.0)*(energy_values[q] -
                              0.5*std::pow(momentum_values[q],2)/density_values[q]);
         if(pressure < eps)
         {
            double drho = density_values[q] - density_average[c];
            double dm = momentum_values[q] - momentum_average[c];
            double dE = energy_values[q] - energy_average[c];
            double a1 = 2.0*drho*dE - dm*dm;
            double b1 = 2.0*drho*(energy_average[c] - eps/(gas_gamma-1.0))
                        + 2.0*density_average[c]*dE
                        - 2.0*momentum_average[c]*dm;
            double c1 = 2.0*density_average[c]*energy_average[c]
                        - momentum_average[c]*momentum_average[c]
                        - 2.0*eps*density_average[c]/(gas_gamma-1.0);
            double D = std::sqrt( std::fabs(b1*b1 - 4.0*a1*c1) );
            double t1 = 0.5*(-b1 - D)/a1;
            double t2 = 0.5*(-b1 + D)/a1;
            double t;
            if(t1 > -1.0e-12 && t1 < 1.0 + 1.0e-12)
               t = t1;
            else if(t2 > -1.0e-12 && t2 < 1.0 + 1.0e-12)
                  t = t2;
            else
            {
               std::cout << "Problem t1, t2 = " << t1 << " " << t2 << "\n";
               std::cout << "eps, rho_min = " << eps << " " << rho_min << "\n";
               std::cout << "theta1 = " << theta1 << "\n";
               std::cout << "pressure = " << pressure << "\n";
               exit(0);
            }
            // t should strictly lie in [0,1]
            t = std::min(1.0, t);
            t = std::max(0.0, t);
            // Need t < 1.0. If t==1 upto machine precision
            // then we are suffering from round off error.
            // In this case we take the cell average value, t=0.
            if(std::fabs(1.0-t) < 1.0e-14) t = 0.0;
            theta2 = std::min(theta2, t);
         }
      }
      
      for(unsigned int i=1; i<dofs_per_cell; ++i)
      {
         density(local_dof_indices[i])  *= theta2;
         momentum(local_dof_indices[i]) *= theta2;
         energy(local_dof_indices[i])   *= theta2;
      }
   }
}


//------------------------------------------------------------------------------
// Compute time step from cfl condition
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cells; ++i)
   {
      double velocity = momentum_average[i] / density_average[i];
      double pressure = (gas_gamma-1.0) * ( energy_average[i] -
                        0.5 * momentum_average[i] * velocity );
      double sonic = std::sqrt ( gas_gamma * pressure / density_average[i] );
      double speed = std::fabs(velocity) + sonic;
      dt = std::min (dt, dx/speed);
   }
   
   dt *= cfl;
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::update (const unsigned int rk_stage)
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

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::output_results () const
{
   // counter to set file name
   static unsigned int c = 0;
   
   DataOut<dim> data_out;
   
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (density, "density");
   
   if(fe.degree <= 1)
      data_out.build_patches (1);
   else
      data_out.build_patches (fe.degree+1);
   
   std::string filename = "sol_" + Utilities::int_to_string(c) + ".gpl";
   std::ofstream output (filename);
   data_out.write_gnuplot (output);

   // save cell average solution
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   std::ofstream fo;
   filename = "avg.gpl";
   fo.open (filename);

   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      Point<dim> x = cell->center();
      double velocity = momentum_average[c] / density_average[c];
      double pressure = (gas_gamma-1.0) * ( energy_average[c] -
                        0.5 * momentum_average[c] * velocity );
      int ind = (is_troubled[c]) ? 1 : 0;
      fo << x(0) << " " 
         << density_average[c] << "  " 
         << velocity << "  " 
         << pressure << "  " 
         << ind << "  "
         << viscosity(c)
         << std::endl;
   }

   fo.close ();

   // increment filename counter
   ++c;
}
//------------------------------------------------------------------------------
// Compute error in solution
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_errors(double& L2_error,
                                       double& H1_error,
                                       double& Linf_error) const
{
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      density,
                                      ExactSolution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree+2),
                                      VectorTools::L2_norm);
   L2_error = difference_per_cell.l2_norm();
   
   VectorTools::integrate_difference (dof_handler,
                                      density,
                                      ExactSolution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree+2),
                                      VectorTools::H1_seminorm);
   H1_error = difference_per_cell.l2_norm();
   
   VectorTools::integrate_difference (dof_handler,
                                      density,
                                      ExactSolution<dim>(),
                                      difference_per_cell,
                                      QGaussLobatto<dim>(fe.degree+2),
                                      VectorTools::Linfty_norm);
   Linf_error = difference_per_cell.l2_norm();
}
//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::run (double& h,
                             int& ndof,
                             double& L2_error,
                             double& H1_error,
                             double& Linf_error)
{
    std::cout << "\n Solving 1-D Euler problem ...\n";

    make_grid_and_dofs();
    assemble_mass_matrix ();
    initialize ();
    compute_averages ();
    identify_troubled_cells ();
    apply_limiter ();
    if(lim_pos) apply_positivity_limiter ();
    output_results ();

    double time = 0.0;
    unsigned int iter = 0;

    std::cout << "Starting the time stepping ... \n";

    while (time < final_time && iter < max_iter)
    {
       density_old  = density;
       momentum_old = momentum;
       energy_old   = energy;
       
       compute_dt ();
       if(time+dt > final_time) dt = final_time - time;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
          compute_viscosity ();
          assemble_rhs ();
          update (rk);
          compute_averages ();
          identify_troubled_cells ();
          apply_limiter ();
          if(lim_pos) apply_positivity_limiter ();
       }
       
       if(iter==0)
       {
          std::cout << "Initial residual = " << residual[0] << " "
                    << residual[1] << " "
                    << residual[2] << std::endl;
          for(unsigned int i=0; i<3; ++i)
             residual0[i] = residual[i];
       }
       
       for(unsigned int i=0; i<3; ++i)
          residual[i] /= residual0[i];
       
      time += dt;
      ++iter;
      if(iter % save_freq == 0) output_results ();
       
      std::cout << "Iter = " << iter << " time = " << time 
                << " Res =" << residual[0] << " " << residual[1] << " "
                << residual[2] << std::endl;
    }
    output_results ();
   
   if(test_case == "smooth") compute_errors (L2_error, H1_error, Linf_error);
   h = dx;
   ndof = dof_handler.n_dofs();
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("degree","0", Patterns::Integer(0,6),
                     "Polynomial degree");
   prm.declare_entry("ncells","100", Patterns::Integer(10,100000),
                     "Number of elements");
   prm.declare_entry("save frequency","100000", Patterns::Integer(0,100000),
                     "How often to save solution");
   prm.declare_entry("test case","sod", 
                     Patterns::Selection("sod|lowd|blast|blastwc|lax|shuosher|sedov|smooth"),
                     "Test case");
   prm.declare_entry("limiter","TVB", 
                     Patterns::Selection("None|TVB|BDF|BSB|visc"),
                     "limiter");
   prm.declare_entry("viscosity","constant",
                     Patterns::Selection("constant|persson"),
                     "artificial viscosity");
   prm.declare_entry("indicator","None",
                     Patterns::Selection("None|density|energy|entropy"),
                     "Shock indicator");
   prm.declare_entry("flux","kfvs", 
                     Patterns::Selection("kfvs|lxf|roe"),
                     "limiter");
   prm.declare_entry("characteristic limiter", "false",
                     Patterns::Bool(), "Characteristic limiter");
   prm.declare_entry("positivity limiter", "false",
                     Patterns::Bool(), "positivity limiter");
   prm.declare_entry("cfl", "1.0",
                     Patterns::Double(0,1.0), "cfl number");
   prm.declare_entry("cip", "0.0",
                     Patterns::Double(0,100.0), "IP constant");
   prm.declare_entry("M", "0.0",
                     Patterns::Double(0,1.0e20), "TVB constant");
   prm.declare_entry("refine","0", Patterns::Integer(0,10),
                     "Number of mesh refinements");
   prm.declare_entry("max iter","1000000000", Patterns::Integer(0,1000000000),
                     "maximum iterations");
}
//------------------------------------------------------------------------------
// Compute convergence rates
//------------------------------------------------------------------------------
void compute_rate(std::vector<double>& h, std::vector<int>& ndof,
                  std::vector<double>& L2_error, std::vector<double>& H1_error,
                  std::vector<double>& Linf_error)
{
   ConvergenceTable   convergence_table;
   unsigned int nrefine = h.size() - 1;
   for(unsigned int i=0; i<=nrefine; ++i)
   {
      convergence_table.add_value("cycle", i);
      convergence_table.add_value("h", h[i]);
      convergence_table.add_value("dofs", ndof[i]);
      convergence_table.add_value("L2", L2_error[i]);
      convergence_table.add_value("H1", H1_error[i]);
      convergence_table.add_value("Linf", Linf_error[i]);
   }

   convergence_table.set_precision("L2", 3);
   convergence_table.set_precision("H1", 3);
   convergence_table.set_precision("Linf", 3);

   convergence_table.set_scientific("h", true);
   convergence_table.set_scientific("L2", true);
   convergence_table.set_scientific("H1", true);
   convergence_table.set_scientific("Linf", true);

   convergence_table.set_tex_caption("h", "$h$");
   convergence_table.set_tex_caption("dofs", "\\# dofs");
   convergence_table.set_tex_caption("L2", "$L^2$-error");
   convergence_table.set_tex_caption("H1", "$H^1$-error");
   convergence_table.set_tex_caption("Linf", "$L_\\infty$-error");

   convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
   convergence_table.evaluate_convergence_rates("H1", ConvergenceTable::reduction_rate_log2);
   convergence_table.evaluate_convergence_rates("Linf", ConvergenceTable::reduction_rate_log2);

   std::cout << std::endl;
   convergence_table.write_text (std::cout);

   std::ofstream error_table_file("error.tex");
   convergence_table.write_tex (error_table_file);

}
//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main (int argc, char** argv)
{
   deallog.depth_console (0);
   {
      ParameterHandler prm;
      declare_parameters (prm);
      if(argc < 2)
      {
         std::cout << "Specify input parameter file\n";
         std::cout << "It should contain following parameters.\n\n";
         prm.print_parameters(std::cout, ParameterHandler::Text);
         return 0;
      }
      bool status = prm.read_input (argv[1], true);
      AssertThrow( status, ExcFileNotOpen(argv[1]) );
      prm.print_parameters(std::cout, ParameterHandler::Text);
      unsigned int degree = prm.get_integer("degree");
      unsigned int nrefine = prm.get_integer("refine");
      std::vector<double> h(nrefine+1), L2_error(nrefine+1), H1_error(nrefine+1),
                          Linf_error(nrefine+1);
      std::vector<int> ndof(nrefine+1);
      for(unsigned int i=0; i<=nrefine; ++i)
      {
         EulerProblem<1> euler_problem(degree, prm);
         euler_problem.run (h[i], ndof[i], L2_error[i], H1_error[i], Linf_error[i]);
         const long int ncells = 2 * prm.get_integer("ncells");
         prm.set("ncells", ncells);
      }
      if(nrefine > 0) compute_rate(h, ndof, L2_error, H1_error, Linf_error);
   }

   return 0;
}
