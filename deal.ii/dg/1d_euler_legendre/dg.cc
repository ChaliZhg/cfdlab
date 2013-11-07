/* 1d DG code for euler equations. This is not particularly efficient code.
   TODO : Use MeshWorker to assemble rhs.
   * Legendre basis functions
   * TVD/TVB limiter for degree = 1,2
   * Characteristic based limiter
   * Positivity preserving limiter
   * Numerical fluxes: Lax-Friedrich, KFVS
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

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {lxf, kfvs};
enum TestCase {sod, blast, blastwc, lax, shuosher, lowd};

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
// Initial condition
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>() {}
   
   virtual void vector_value (const Point<dim>   &p,
                              Vector<double>& values) const;
   string test_case;
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
         values(0) = d_left; // left density
         values(1) = u_left; // left velocity
         values(2) = p_left; // left pressure
      }
      else
      {
         values(0) = d_right; // right density
         values(1) = u_right; // right velocity
         values(2) = p_right; // right pressure
      }
   }
   else if(test_case == "blast")
   {
      if(p[0] < 0.1)
      {
         values(0) = 1; // left density
         values(1) = 0; // left velocity
         values(2) = 1000; // left pressure
      }
      else if(p[0] > 0.9)
      {
         values(0) = 1; // right density
         values(1) = 0; // right velocity
         values(2) = 100; // right pressure
      }
      else
      {
         values(0) = 1; // right density
         values(1) = 0; // right velocity
         values(2) = 0.01; // right pressure
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
   void run ();
   
private:
   void make_grid_and_dofs ();
   void initialize ();
   void assemble_mass_matrix ();
   void assemble_rhs ();
   void compute_averages ();
   void compute_dt ();
   void apply_limiter ();
   void apply_positivity_limiter ();
   void update (const unsigned int rk_stage);
   void output_results () const;
   
   unsigned int         n_cells;
   string               test_case;
   double               dt;
   double               dx;
   double               cfl;
   double               final_time;
   double               min_residue;
   unsigned int         n_rk_stages;
   FluxType             flux_type;
   bool                 lim_char, lim_pos;
   bool                 lbc_reflect, rbc_reflect;
   unsigned int         save_freq;
   
   
   Triangulation<dim>   triangulation;
   FE_DGP<dim>          fe;
   DoFHandler<dim>      dof_handler;
   
   SparsityPattern      sparsity_pattern;
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
   
   std::vector< Vector<double> > density_average;
   std::vector< Vector<double> > momentum_average;
   std::vector< Vector<double> > energy_average;
   
   std::vector<double> residual;
   std::vector<double> residual0;
   
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
   double M = prm.get_double("M");
   save_freq= prm.get_integer("save frequency");
   
   n_rk_stages = 3;
   flux_type = kfvs;
   
   lbc_reflect = rbc_reflect = false;
   
   if(test_case == "sod")
   {
      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.5;
      final_time = 0.2;
      min_residue= 1.0e20;
      
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
      min_residue= 1.0e20;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      lbc_reflect = rbc_reflect = true;
   }
   else if(test_case == "shuosher")
   {
      xmin    = -5;
      xmax    = 5;
      final_time = 1.8;
      min_residue= 1.0e20;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
   }
   else if(test_case == "lax")
   {
      xmin    = -5;
      xmax    = 5;
      xmid    = 0.0;
      final_time = 1.3;
      min_residue= 1.0e20;
      
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
      min_residue= 1.0e20;
      
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
      min_residue= 1.0e20;
      cfl      = 0.8;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      
      d_left  = 1.0;
      d_right = 1.0;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1000.0;
      p_right = 0.01;
      save_freq = 10;
      M = 500.0;
      
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

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);
   
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
   
    density_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
    momentum_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
    energy_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
   
    residual.resize(3, 1.0);
    residual0.resize(3);
   
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
// Flux for NS equation
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
// KFVS flux for navier-stokes
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
// Compute flux across cell faces
//------------------------------------------------------------------------------
void numerical_flux (const FluxType& flux_type,
                     Vector<double>& left_state,
                     Vector<double>& right_state,
                     Vector<double>& flux)
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
   
   // for getting neighbor cell solution using trapezoidal rule
   std::vector<double>  density_values_n  (2);
   std::vector<double>  momentum_values_n (2);
   std::vector<double>  energy_values_n   (2);

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

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            euler_flux(density_values[q_point], momentum_values[q_point],
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
       
       // Computation of flux at cell boundaries
       Vector<double> lf_left_state(3), lf_right_state(3);

        // left face flux
        // right state is from current cell
       lf_right_state(0) = density_values [0];
       lf_right_state(1) = momentum_values[0];
       lf_right_state(2) = energy_values  [0];
       
       if(c==0)
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
//             lf_left_state(0) = d_left;
//             lf_left_state(1) = d_left * u_left;
//             lf_left_state(2) = p_left/(gas_gamma-1.0) + 0.5 * d_left * std::pow(u_left,2);
          }
       }
       else
       {
          // get left cell dof indices
          fe_values_neighbor.reinit (cell->neighbor(0));
          
          fe_values_neighbor.get_function_values (density,  density_values_n);
          fe_values_neighbor.get_function_values (momentum, momentum_values_n);
          fe_values_neighbor.get_function_values (energy,   energy_values_n);
          
          lf_left_state(0) = density_values_n [1];
          lf_left_state(1) = momentum_values_n[1];
          lf_left_state(2) = energy_values_n  [1];
          
       }
       
       Vector<double> left_flux(3);
       numerical_flux (flux_type, lf_left_state, lf_right_state, left_flux);
       
       std::vector<double> lf_left_prim = con2prim(lf_left_state);
       std::vector<double> lf_right_prim = con2prim(lf_right_state);
       
       // right face flux
       Vector<double> rf_left_state(3), rf_right_state(3);
       // left state is from current cell
       rf_left_state(0) = density_values [n_q_points-1];
       rf_left_state(1) = momentum_values[n_q_points-1];
       rf_left_state(2) = energy_values  [n_q_points-1];
       
       if(c==triangulation.n_cells()-1)
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
//             rf_right_state(0) = d_right;
//             rf_right_state(1) = d_right * u_right;
//             rf_right_state(2) = p_right/(gas_gamma-1.0) +
//                                 0.5 * d_right * std::pow(u_right,2);
          }
       }
       else
       {          
          // get right cell to right face
          fe_values_neighbor.reinit (cell->neighbor(1));
          
          fe_values_neighbor.get_function_values (density,  density_values_n);
          fe_values_neighbor.get_function_values (momentum, momentum_values_n);
          fe_values_neighbor.get_function_values (energy,   energy_values_n);
          
          rf_right_state(0) = density_values_n [0];
          rf_right_state(1) = momentum_values_n[0];
          rf_right_state(2) = energy_values_n  [0];
       }
       
       Vector<double> right_flux(3);
       numerical_flux (flux_type, rf_left_state, rf_right_state, right_flux);
       
       std::vector<double> rf_left_prim = con2prim(rf_left_state);
       std::vector<double> rf_right_prim = con2prim(rf_right_state);
       
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
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_averages ()
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
      
      density_average[c] = 0.0;
      momentum_average[c] = 0.0;
      energy_average[c] = 0.0;
      for(unsigned int point=0; point<n_q_points; ++point)
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            density_average[c](0) += density(local_dof_indices[i]) * 
                                     fe_values.shape_value (i, point) *
                                     fe_values.JxW (point);
            momentum_average[c](0) += momentum(local_dof_indices[i]) * 
                                      fe_values.shape_value (i, point) *
                                      fe_values.JxW (point);
            energy_average[c](0) += energy(local_dof_indices[i]) * 
                                    fe_values.shape_value (i, point) *
                                    fe_values.JxW (point);
            
            if(fe.degree >= 1) // compute average gradient
            {
               density_average[c](1) += density(local_dof_indices[i]) * 
                                        fe_values.shape_grad (i, point)[0] *
                                        fe_values.JxW (point);
               momentum_average[c](1) += momentum(local_dof_indices[i]) * 
                                         fe_values.shape_grad (i, point)[0] *
                                         fe_values.JxW (point);
               energy_average[c](1) += energy(local_dof_indices[i]) * 
                                       fe_values.shape_grad (i, point)[0] *
                                       fe_values.JxW (point);
            }
         }
      
      density_average[c]  /= dx;
      momentum_average[c] /= dx;
      energy_average[c]   /= dx;
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
// currently works only for P1/P2 case
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_limiter ()
{
   if(fe.degree == 0) return;
   
   Assert (fe.degree<3, ExcIndexRange(fe.degree, 1, 3));
   
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
   {
      fe_values.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      fe_values.get_function_values(density, density_face_values);
      fe_values.get_function_values(momentum, momentum_face_values);
      fe_values.get_function_values(energy, energy_face_values);
      
      if(c==0)
      {
         density_left = density_average[c](0);
         if(lbc_reflect)
            momentum_left = -momentum_average[c](0);
         else
            momentum_left = momentum_average[c](0);
         energy_left = energy_average[c](0);
         
         density_right = density_average[c+1](0);
         momentum_right = momentum_average[c+1](0);
         energy_right = energy_average[c+1](0);
      }
      else if(c == n_cells-1)
      {
         density_left = density_average[c-1](0);
         momentum_left = momentum_average[c-1](0);
         energy_left = energy_average[c-1](0);
         
         density_right = density_average[c](0);
         if(rbc_reflect)
            momentum_right = -momentum_average[c](0);
         else
            momentum_right = momentum_average[c](0);
         energy_right = energy_average[c](0);
      }
      else
      {
         density_left = density_average[c-1](0);
         momentum_left = momentum_average[c-1](0);
         energy_left = energy_average[c-1](0);
         
         density_right = density_average[c+1](0);
         momentum_right = momentum_average[c+1](0);
         energy_right = energy_average[c+1](0);
      }
      
      // density
      db[0] = density_average[c](0) - density_left;
      df[0] = density_right - density_average[c](0);
      DB[0] = density_average[c](0) - density_face_values[0];
      DF[0] = density_face_values[1] - density_average[c](0);
      
      // momentum
      db[1] = momentum_average[c](0) - momentum_left;
      df[1] = momentum_right - momentum_average[c](0);
      DB[1] = momentum_average[c](0) - momentum_face_values[0];
      DF[1] = momentum_face_values[1] - momentum_average[c](0);
      
      // energy
      db[2] = energy_average[c](0) - energy_left;
      df[2] = energy_right - energy_average[c](0);
      DB[2] = energy_average[c](0) - energy_face_values[0];
      DF[2] = energy_face_values[1] - energy_average[c](0);

      double R[n_var][n_var], Ri[n_var][n_var];
      if(lim_char)
      {
         EigMat(density_average[c](0), 
                momentum_average[c](0), 
                energy_average[c](0), R, Ri);
         Multi(Ri, db);
         Multi(Ri, df);
         Multi(Ri, DB);
         Multi(Ri, DF);
      }

      for(unsigned int i=0; i<n_var; ++i)
      {
         dl[i] = minmod (DB[i], db[i], df[i]);
         dr[i] = minmod (DF[i], db[i], df[i]);
      }
      
      if(fe.degree==1) // dl = dr
      {
         if(lim_char) Multi(R, dr);
         density(local_dof_indices[1])  = dr[0] / fe_values.shape_value(1,1);
         momentum(local_dof_indices[1]) = dr[1] / fe_values.shape_value(1,1);
         energy(local_dof_indices[1])   = dr[2] / fe_values.shape_value(1,1);
      }
      else if(fe.degree==2)
      {
         if(lim_char)
         {
            Multi(R, dl);
            Multi(R, dr);
         }
         double phi1_l = fe_values.shape_value(1,0);
         double phi1_r = fe_values.shape_value(1,1);
         double phi2_l = fe_values.shape_value(2,0);
         double phi2_r = fe_values.shape_value(2,1);
         double det = - phi1_l * phi2_r + phi1_r * phi2_l;
         density(local_dof_indices[1])  = (phi2_r * dl[0] + phi2_l * dr[0])/det;
         momentum(local_dof_indices[1]) = (phi2_r * dl[1] + phi2_l * dr[1])/det;
         energy(local_dof_indices[1])   = (phi2_r * dl[2] + phi2_l * dr[2])/det;
         
         density(local_dof_indices[2])  = (-phi1_r * dl[0] - phi1_l * dr[0])/det;
         momentum(local_dof_indices[2]) = (-phi1_r * dl[1] - phi1_l * dr[1])/det;
         energy(local_dof_indices[2])   = (-phi1_r * dl[2] - phi1_l * dr[2])/det;
      }
      
   }
}

//------------------------------------------------------------------------------
// Apply positivity limiter
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::apply_positivity_limiter ()
{
   if(fe.degree == 0) return;
   
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
      double velocity = momentum_average[c](0) / density_average[c](0);
      double pressure = (gas_gamma-1.0) * ( energy_average[c](0) -
                                           0.5 * momentum_average[c](0) * velocity );
      eps = std::min(eps, density_average[c](0));
      eps = std::min(eps, pressure);
   }
   if(eps < 0.0)
   {
      std::cout << "Negative states\n";
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
      
      double theta1 =
         std::min(std::fabs((density_average[c](0) - eps)/(density_average[c](0) - rho_min)), 1.0);
      
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
            double drho = density_average[c](0) - density_values[q];
            double dm = momentum_average[c](0) - momentum_values[q];
            double dE = energy_average[c](0) - energy_values[q];
            double a1 = 2.0*drho*dE - dm*dm;
            double b1 = 2.0*energy_values[q]*drho + 2.0*density_values[q]*dE
                        - 2.0*momentum_values[q]*dm - 2.0*eps*drho/(gas_gamma-1.0);
            double c1 = 2.0*density_values[q]*energy_values[q]
                        - std::pow(momentum_values[q],2.0)
                        - 2.0*eps*density_values[q]/(gas_gamma-1.0);
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
            theta2 = std::min(theta2, 1.0-t);
               std::cout << "Problem t1, t2 = " << t1 << " " << t2 << "\n";
               std::cout << "eps, rho_min = " << eps << " " << rho_min << "\n";
               std::cout << "theta1,theta2 = " << theta1 << " " << theta2 << "\n";
               std::cout << "pressure = " << pressure << "\n";
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
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cells; ++i)
   {
      double velocity = momentum_average[i](0) / density_average[i](0);
      double pressure = (gas_gamma-1.0) * ( energy_average[i](0) -
                        0.5 * momentum_average[i](0) * velocity );
      double sonic = std::sqrt ( gas_gamma * pressure / density_average[i](0) );
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
   static unsigned int c = 0;
   
   Vector<double> velocity(dof_handler.n_dofs());
   Vector<double> pressure(dof_handler.n_dofs());
   
   // Compute velocity and pressure
   for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   {
      velocity(i) = momentum(i) / density(i);
      pressure(i) = (gas_gamma-1.0) * (energy(i) -
                                       0.5 * momentum(i) * velocity(i));
   }
   
   DataOut<dim> data_out;
   
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (density, "density");
   data_out.add_data_vector (velocity, "velocity");
   data_out.add_data_vector (pressure, "pressure");
   
   if(fe.degree <= 1)
      data_out.build_patches (1);
   else
      data_out.build_patches (5);
   
   std::string filename = "sol_" + Utilities::int_to_string(c) + ".gpl";
   std::ofstream output (filename);
   data_out.write_gnuplot (output);
   ++c;
   
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void EulerProblem<dim>::run ()
{
    std::cout << "\n Solving 1-D NS problem ...\n";

    make_grid_and_dofs();
    assemble_mass_matrix ();
    initialize ();
    output_results ();
    compute_averages ();
    apply_limiter ();

    double time = 0.0;
    unsigned int iter = 0;

    while (time < final_time || ( residual[0] > min_residue &&
           residual[1] > min_residue && residual[2] > min_residue))
    {
       density_old  = density;
       momentum_old = momentum;
       energy_old   = energy;
       
       compute_dt ();
       if(time+dt > final_time) dt = final_time - time;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
         assemble_rhs ();
         update (rk);
         compute_averages ();
         apply_limiter ();
         if(lim_pos) apply_positivity_limiter ();
       }
       
       if(iter==0)
       {
          std::cout << "Initial residual = " << residual[0] << " "
                    << residual[1] << " "
                    << residual[2] << endl;
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
                << residual[2] << endl;
    }
    output_results ();
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("degree","0", Patterns::Integer(0,2),
                     "Polynomial degree");
   prm.declare_entry("ncells","100", Patterns::Integer(10,100000),
                     "Number of elements");
   prm.declare_entry("save frequency","100000", Patterns::Integer(0,100000),
                     "How often to save solution");
   prm.declare_entry("test case","sod", 
                     Patterns::Selection("sod|lowd|blast|blastwc|lax|shuosher"),
                     "Test case");
   prm.declare_entry("characteristic limiter", "false",
                     Patterns::Bool(), "Characteristic limiter");
   prm.declare_entry("positivity limiter", "false",
                     Patterns::Bool(), "positivity limiter");
   prm.declare_entry("cfl", "1.0",
                     Patterns::Double(0,1.0), "cfl number");
   prm.declare_entry("M", "0.0",
                     Patterns::Double(0,1.0e20), "TVB constant");
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main (int argc, char** argv)
{
   deallog.depth_console (0);
   {
      if(argc < 2)
      {
         std::cout << "Specify input parameter file\n";
         return 0;
      }
      ParameterHandler prm;
      declare_parameters (prm);
      prm.read_input (argv[1]);
      prm.print_parameters(std::cout, ParameterHandler::Text);
      unsigned int degree = prm.get_integer("degree");
      EulerProblem<1> euler_problem(degree, prm);
      euler_problem.run ();
   }

    return 0;
}
