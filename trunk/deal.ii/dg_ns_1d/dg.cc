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
#include <numerics/fe_field_function.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

// Number of variables: mass, momentum and energy
const unsigned int n_var = 3;
double gas_gamma;
double gas_const;
double mu_ref;
double T_ref;
double Pr;
double               d_left, u_left, p_left;
double               d_right, u_right, p_right;
double               xmin, xmax, xmid;

const double NIPG = -1.0;
const double SIPG =  1.0;

// Factor in Maxwell distribution
double alpha;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {lxf, kfvs};
enum TestCase {sod, shock_structure};
enum ViscScheme {none, split, avg};

//------------------------------------------------------------------------------
// Viscosity coefficient as function of temperature
//------------------------------------------------------------------------------
double viscosity (const double& T)
{
   return mu_ref * pow( T/T_ref, 0.8);
}

double viscosity (const Vector<double>& state)
{
   double velocity = state(1) / state(0);
   double pressure = (gas_gamma - 1.0) * (state(2) - 0.5 * state(1) * velocity);
   
   double T   = pressure / state(0) / gas_const;
   double mu  = viscosity (T);
   return mu;
}
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
   
   virtual void vector_value (const Point<dim>   &p,
                              Vector<double>& values) const;
};

// Initial condition for density, velocity, pressure
template<int dim>
void InitialCondition<dim>::vector_value (const Point<dim>   &p,
                                          Vector<double>& values) const
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

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class NSProblem
{
public:
    NSProblem (unsigned int degree, TestCase test_case);
    void run ();

private:
    void make_grid_and_dofs ();
    void initialize ();
    void assemble_mass_matrix ();
    void assemble_rhs ();
    void compute_averages ();
    void compute_dt ();
    void apply_limiter ();
    void update (const unsigned int rk_stage);
    void output_results () const;
   
    TestCase             test_case;
    unsigned int         n_cells;
    double               dt;
    double               dx;
    double               cfl;
    double               final_time;
    double               min_residue;
    unsigned int         n_rk_stages;
    FluxType             flux_type;
    ViscScheme           visc_scheme;
    double               ip_scheme;
    double               CIP;


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
   
    std::vector< Vector<double> > density_average;
    std::vector< Vector<double> > momentum_average;
    std::vector< Vector<double> > energy_average;
   
    std::vector<double> residual;
    std::vector<double> residual0;

};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
NSProblem<dim>::NSProblem (unsigned int degree, TestCase test_case) :
    test_case (test_case),
    fe (degree),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
   

   n_cells = 200;

   n_rk_stages = 3;
   flux_type = kfvs;
   ip_scheme = SIPG;
   CIP = 10.0;
   
   if(test_case == sod)
   {
      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.5;
      final_time = 0.2;
      min_residue= 1.0e20;
      cfl     = 0.1;
      
      gas_gamma = 1.4;
      gas_const = 1.0;
      Pr = 1.0;
      
      d_left  = 1.0;
      d_right = 0.125;
      
      u_left  = 0.0;
      u_right = 0.0;
      
      p_left  = 1.0;
      p_right = 0.1;
      
      mu_ref = 0.0;
      T_ref  = p_left / d_left / gas_const; // not used
      
      visc_scheme = none;
   }
   else if (test_case == shock_structure)
   {
      xmin    = -0.25;
      xmax    =  0.0;
      xmid    = 0.5 * ( xmin + xmax );
      final_time = 0.5;
      min_residue= 1.0e-8;
      cfl = 0.025;
      
      gas_gamma = 5.0/3.0;
      gas_const = 0.5;
      Pr = 2.0/3.0;
      
      double mach_left = 1.5;
      double M2 = pow(mach_left, 2);
      
      d_left = 1.0;
      u_left = 1.0;
      p_left = 1.0/gas_gamma/M2;
      
      d_right= (gas_gamma+1.0)*M2/(2.0+(gas_gamma-1.0)*M2)*d_left;
      u_right= ((gas_gamma-1.0)/(gas_gamma+1.0)+2.0/(gas_gamma+1.0)/M2)*u_left;
      p_right= (2.0*gas_gamma/(gas_gamma+1.0)*M2-(gas_gamma-1)/(gas_gamma+1.0))*p_left;
      
      mu_ref= 0.0005;
      T_ref = p_left / d_left / gas_const;
      
      visc_scheme = split;
   }
   else
   {
      std::cout << "Unknown test case\n";
   }
   
   dx      = (xmax - xmin) / n_cells;

   double d = 1.0/( 1.0/(gas_gamma - 1.0) - dim/2.0 );
   alpha = std::sqrt(2.0*M_PI) * gamma(1.0 + 1.0/d);

}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
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
    face_flux.resize(triangulation.n_active_cells()+1, Vector<double>(n_var));
   
    density_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
    momentum_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
    energy_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
   
    residual.resize(3, 1.0);
    residual0.resize(3);
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
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
   Vector<double> initial_value(n_var);
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

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// compute tau and q
// state = conserved
// grad  = conserved grad
//------------------------------------------------------------------------------
void compute_nsterms (const Vector<double>& state,
                      const Vector<double>& grad,
                      double& tau,
                      double& q)
{
   double velocity = state(1) / state(0);
   double pressure = (gas_gamma - 1.0) * (state(2) - 0.5 * state(1) * velocity);
   
   double T   = pressure / state(0) / gas_const;
   double mu  = viscosity (T);
   double u_x = (grad(1) - velocity * grad(0)) / state(0);
   tau = 4.0 * mu * u_x / 3.0;
   
   double T_x = ((gas_gamma-1.0)/gas_const * 
                 (grad(2) - 0.5 * velocity * grad(1) - 0.5 * state(1) * u_x)
                 - T * grad(0)) / state(0);
   double k   = mu * gas_gamma * gas_const / (gas_gamma - 1.0) / Pr;
   q   = -k * T_x;
}
//------------------------------------------------------------------------------
// Flux for NS equation
//------------------------------------------------------------------------------
void ns_flux (const double& density,
              const double& momentum,
              const double& energy,
              const Tensor<1,1>& density_grad,
              const Tensor<1,1>& momentum_grad,
              const Tensor<1,1>& energy_grad,
              Vector<double>& flux)
{   
   double velocity = momentum / density;
   double pressure = (gas_gamma - 1.0) * (energy - 0.5 * momentum * velocity);
   flux(0) = momentum;
   flux(1) = pressure + momentum * velocity;
   flux(2) = (energy + pressure) * velocity;
   
   Vector<double> state (n_var);
   state(0) = density;
   state(1) = momentum;
   state(2) = energy;
   Vector<double> grad (n_var);
   grad(0) = density_grad[0];
   grad(1) = momentum_grad[0];
   grad(2) = energy_grad[0];
   
   double tau, q;
   compute_nsterms (state, grad, tau, q);
      
   flux(1) -= tau;
   flux(2) += -tau * velocity + q;
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
void kfvs_split_flux (unsigned int mode,
                      const std::vector<double>& prim,
                      const double& tau,
                      const double& q,
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
   
   if(mode==0)
   {
      // inviscid flux
      flux[0] = prim[0] * fact;
      flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A +
                prim[0] * prim[1] * B;
      flux[2] = prim[1] * (E + prim[2]) * A +
                (E + 0.5 * prim[2]) * B;
   }
   else
      flux[0] = flux[1] = flux[2] = 0.0;
   
   static const double g1 = (gas_gamma - 1.0) / gas_gamma;
   static const double g2 = (3.0 * gas_gamma - 1.0) / (4.0 * (gas_gamma - 1.0) );
   static const double g3 = (3.0 - gas_gamma) / (2.0 * gas_gamma);
   
   // viscous flux
   flux[0] += prim[0] * B / prim[2] * ( -0.5 * tau -
                                       g1 * prim[1] * beta * q );
   flux[1] += -tau * A + g1 * prim[0] * B * q / prim[2];
   flux[2] += (-prim[1] * tau + q) * A +
              B * (-g2 * tau - g3 * prim[1] * beta * q);
}
//------------------------------------------------------------------------------
// KFVS flux for navier-stokes
//------------------------------------------------------------------------------
void KFVSFlux (const Vector<double>& left_state,
               const Vector<double>& right_state,
               const Vector<double>& left_grad,
               const Vector<double>& right_grad,
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
   
   double tau_l, q_l;
   compute_nsterms (left_state, left_grad, tau_l, q_l);

   double tau_r, q_r;
   compute_nsterms (right_state, right_grad, tau_r, q_r);
   
   std::vector<double> flux_pos (n_var);
   std::vector<double> flux_neg (n_var);
   
   kfvs_split_flux (0, left,  tau_l, q_l, +1, flux_pos);
   kfvs_split_flux (0, right, tau_r, q_r, -1, flux_neg);
   
   for(unsigned int i=0; i<n_var; ++i)
      flux(i) = flux_pos[i] + flux_neg[i];
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void numerical_flux (const FluxType& flux_type,
                     Vector<double>& left_state,
                     Vector<double>& right_state,
                     Vector<double>& left_grad,
                     Vector<double>& right_grad,
                     Vector<double>& flux)
{
   switch (flux_type) 
   {
      case lxf:
         LaxFlux (left_state, right_state, flux);
         break;
         
      case kfvs:
         KFVSFlux (left_state, right_state, left_grad, right_grad, flux);
         break;
         
      default:
         cout << "Unknown flux_type !!!\n";
         abort ();
   }
}
//------------------------------------------------------------------------------
//Compute entropy variable given conserved state
//------------------------------------------------------------------------------
Vector<double> entropy_var (const Vector<double>& state)
{
   Vector<double> evar(n_var);
   
   double velocity = state(1) / state(0);
   double pressure = (gas_gamma-1.0) * (state(2) - 
                                        0.5 * state(1) * velocity );
   double T = pressure / (gas_const * state(0));
   
   double ss = state(0) / std::pow(gas_const*T, 1.0/(gas_gamma-1)) / alpha;
   evar(0) = std::log(ss) - 0.5 * velocity * velocity / (gas_const * T);
   evar(1) = velocity / (gas_const * T);
   evar(2) = -1.0 / (gas_const * T);
   
   return evar;
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
// convert entropy gradient to conserved gradient
// grad(U) = tilde{A}_0 grad(V)
//------------------------------------------------------------------------------
void convert_con (const std::vector<double>& prim, Vector<double>& grad)
{
   Vector<double> vec = grad;
   FullMatrix<double> mat(n_var, n_var);
   
   double E = prim[2]/(gas_gamma-1) + 0.5 * prim[0] * std::pow(prim[1],2);
   
   mat(0,0) = prim[0];
   mat(0,1) = prim[0] * prim[1];
   mat(0,2) = E;
   
   mat(1,1) = prim[2] + prim[0] * prim[1] * prim[1];
   mat(1,2) = (E + prim[2]) * prim[1];
   
   mat(2,2) = gas_gamma * prim[2] * (E + 0.5 * prim[0] * prim[1] * prim[1]) /
              (prim[0] * (gas_gamma - 1)) +
               0.25 * prim[0] * std::pow(prim[1], 4);
   
   //symmetry
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);
   
   grad = 0.0;
   for(unsigned int i=0; i<n_var; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         grad(i) += mat(i,j) * vec(j);
}

//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::assemble_rhs ()
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
   
    std::vector< Tensor<1,dim> >  density_grad (n_q_points);
    std::vector< Tensor<1,dim> >  momentum_grad (n_q_points);
    std::vector< Tensor<1,dim> >  energy_grad (n_q_points);
   
   // for getting neighbor cell solution using trapezoidal rule
   std::vector<double>  density_values_n  (2);
   std::vector<double>  momentum_values_n (2);
   std::vector<double>  energy_values_n   (2);
   
   std::vector< Tensor<1,dim> >  density_grad_n (2);
   std::vector< Tensor<1,dim> >  momentum_grad_n (2);
   std::vector< Tensor<1,dim> >  energy_grad_n (2);


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

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            ns_flux(density_values[q_point], momentum_values[q_point], 
                    energy_values[q_point], density_grad[q_point],
                    momentum_grad[q_point], energy_grad[q_point], flux);
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

        // left face flux
        // right state is from current cell
       lf_right_state(0) = density_values [0];
       lf_right_state(1) = momentum_values[0];
       lf_right_state(2) = energy_values  [0];
       lf_right_grad(0)  = density_grad [0][0];
       lf_right_grad(1)  = momentum_grad[0][0];
       lf_right_grad(2)  = energy_grad  [0][0];
       
       if(c==0)
       {
          lf_left_state(0) = d_left;
          lf_left_state(1) = d_left * u_left;
          lf_left_state(2) = p_left/(gas_gamma-1.0) + 0.5 * d_left * std::pow(u_left,2);
          lf_left_grad  = 0.0;
       }else
       {
          // get left cell dof indices
          fe_values_neighbor.reinit (cell->neighbor(0));
          
          fe_values_neighbor.get_function_values (density,  density_values_n);
          fe_values_neighbor.get_function_values (momentum, momentum_values_n);
          fe_values_neighbor.get_function_values (energy,   energy_values_n);
          
          fe_values_neighbor.get_function_gradients (density,  density_grad_n);
          fe_values_neighbor.get_function_gradients (momentum, momentum_grad_n);
          fe_values_neighbor.get_function_gradients (energy,   energy_grad_n);
          
          lf_left_state(0) = density_values_n [1];
          lf_left_state(1) = momentum_values_n[1];
          lf_left_state(2) = energy_values_n  [1];
          
          lf_left_grad(0) = density_grad_n [1][0];
          lf_left_grad(1) = momentum_grad_n[1][0];
          lf_left_grad(2) = energy_grad_n  [1][0];
       }
       
       if(visc_scheme == avg) //average gradients
       {
          Vector<double> avg_grad (n_var);
          for(unsigned int i=0; i<n_var; ++i)
             avg_grad(i) = 0.5*(lf_left_grad(i) + lf_right_grad(i));
          lf_left_grad = avg_grad;
          lf_right_grad = avg_grad;
       }
       
       Vector<double> left_flux(3);
       numerical_flux (flux_type, lf_left_state, lf_right_state, 
                       lf_left_grad, lf_right_grad, left_flux);
       
       Vector<double> lf_left_entropy = entropy_var (lf_left_state);
       Vector<double> lf_right_entropy = entropy_var (lf_right_state);
       Vector<double> lf_entropy_jump(n_var);
       Vector<double> lf_state_jump(n_var);
       for(unsigned int i=0; i<n_var; ++i)
       {
          lf_entropy_jump(i) = lf_right_entropy(i) - lf_left_entropy(i);
          lf_state_jump (i) = lf_right_state (i) - lf_left_state (i);
       }
       std::vector<double> lf_left_prim = con2prim(lf_left_state);
       std::vector<double> lf_right_prim = con2prim(lf_right_state);
       
       double lf_left_nu = viscosity (lf_left_state) / lf_left_state(0);
       double lf_right_nu = viscosity (lf_right_state) / lf_right_state(0);
       double lf_nu = 0.5*(lf_left_nu + lf_right_nu);

       
       // right face flux
       Vector<double> rf_left_state(3), rf_right_state(3);
       Vector<double> rf_left_grad(3), rf_right_grad(3);
       // left state is from current cell
       rf_left_state(0) = density_values [n_q_points-1];
       rf_left_state(1) = momentum_values[n_q_points-1];
       rf_left_state(2) = energy_values  [n_q_points-1];
       rf_left_grad(0)  = density_grad [n_q_points-1][0];
       rf_left_grad(1)  = momentum_grad[n_q_points-1][0];
       rf_left_grad(2)  = energy_grad  [n_q_points-1][0];
       
       if(c==triangulation.n_cells()-1)
       {
          rf_right_state(0) = d_right;
          rf_right_state(1) = d_right * u_right;
          rf_right_state(2) = p_right/(gas_gamma-1.0) + 0.5 * d_right * std::pow(u_right,2);
          rf_right_grad  = 0.0;
       }else
       {          
          // get right cell to right face
          fe_values_neighbor.reinit (cell->neighbor(1));
          
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
       }
       
       if(visc_scheme == avg)
       {
          Vector<double> avg_grad (n_var);
          for(unsigned int i=0; i<n_var; ++i)
             avg_grad(i) = 0.5*(rf_left_grad(i) + rf_right_grad(i));
          rf_left_grad = avg_grad;
          rf_right_grad = avg_grad;
       }
       
       Vector<double> right_flux(3);
       numerical_flux (flux_type, rf_left_state, rf_right_state, 
                       rf_left_grad, rf_right_grad, right_flux);
       
       Vector<double> rf_left_entropy = entropy_var (rf_left_state);
       Vector<double> rf_right_entropy = entropy_var (rf_right_state);
       Vector<double> rf_entropy_jump(n_var);
       Vector<double> rf_state_jump(n_var);
       for(unsigned int i=0; i<n_var; ++i)
       {
          rf_entropy_jump(i) = rf_right_entropy(i) - rf_left_entropy(i);
          rf_state_jump (i) = rf_right_state (i) - rf_left_state (i);
       }
       std::vector<double> rf_left_prim = con2prim(rf_left_state);
       std::vector<double> rf_right_prim = con2prim(rf_right_state);

       double rf_left_nu = viscosity (rf_left_state) / rf_left_state(0);
       double rf_right_nu = viscosity (rf_right_state) / rf_right_state(0);
       double rf_nu = 0.5*(rf_left_nu + rf_right_nu);
       
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
           
           if(visc_scheme != none)
           {
              for(unsigned int j=0; j<n_var; ++j)
              {
                 Vector<double> phi_grad(n_var);
                 double tau, q;
                 
                 //left face
                 phi_grad = 0.0;
                 phi_grad (j) = fe_values.shape_grad(i,0)[0];
                 convert_con (lf_right_prim, phi_grad);
                 compute_nsterms (lf_right_state, phi_grad, tau, q);
                 std::vector<double> flux_p_left(n_var);
                 if(visc_scheme == split)
                    kfvs_split_flux (1, lf_right_prim, tau, q, -1, flux_p_left);
                 else
                 {
                    std::vector<double> fp(n_var), fm(n_var);
                    kfvs_split_flux (1, lf_right_prim, tau, q, -1, fm);
                    kfvs_split_flux (1, lf_left_prim, tau, q, +1, fp);
                    for(unsigned int k=0; k<n_var; ++k)
                       flux_p_left[k] = 0.5*(fm[k] + fp[k]);
                 }
                 
                 double left_stab = 0.0;
                 for(unsigned int k=0; k<n_var; ++k)
                    left_stab += flux_p_left[k] * lf_entropy_jump(k);
                 
                 //right face
                 phi_grad = 0.0;
                 phi_grad (j) = fe_values.shape_grad(i,n_q_points-1)[0];
                 convert_con (rf_left_prim, phi_grad);
                 compute_nsterms (rf_left_state, phi_grad, tau, q);
                 std::vector<double> flux_p_right(n_var);
                 if(visc_scheme == split)
                    kfvs_split_flux (1, rf_left_prim, tau, q, +1, flux_p_right);
                 else
                 {
                    std::vector<double> fp(n_var), fm(n_var);
                    kfvs_split_flux (1, rf_left_prim, tau, q, +1, fp);
                    kfvs_split_flux (1, rf_right_prim, tau, q, -1, fm);
                    for(unsigned int k=0; k<n_var; ++k)
                       flux_p_right[k] = 0.5*(fm[k] + fp[k]);
                 }
                 
                 double right_stab = 0.0;
                 for(unsigned int k=0; k<n_var; ++k)
                    right_stab += flux_p_right[k] * rf_entropy_jump(k);
                 
                 double lf_delta = CIP * lf_nu / dx;
                 double rf_delta = CIP * rf_nu / dx;
                 
                 if(j==0)
                 {
                    cell_rhs_density(i) += ip_scheme * (left_stab + right_stab);
                    cell_rhs_density(i) += -lf_delta * 
                                            lf_state_jump(0) *
                                            fe_values.shape_value(i, 0)
                                            +
                                            rf_delta *
                                            rf_state_jump(0) *
                                            fe_values.shape_value(i,n_q_points-1);
                 }
                 else if(j==1)
                 {
                    cell_rhs_momentum(i) += ip_scheme * (left_stab + right_stab);
                    cell_rhs_momentum(i) += -lf_delta * 
                                             lf_state_jump(1) *
                                             fe_values.shape_value(i, 0)
                                             +
                                             rf_delta *
                                             rf_state_jump(1) *
                                             fe_values.shape_value(i,n_q_points-1);
                 }else
                 {
                    cell_rhs_energy(i) += ip_scheme * (left_stab + right_stab);
                    cell_rhs_energy(i) += -lf_delta * 
                                           lf_state_jump(2) *
                                           fe_values.shape_value(i, 0)
                                           +
                                           rf_delta *
                                           rf_state_jump(2) *
                                           fe_values.shape_value(i,n_q_points-1);
                 }
              }
           }
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
void NSProblem<dim>::compute_averages ()
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
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::apply_limiter ()
{
   Assert (fe.degree==1, ExcIndexRange(fe.degree, 1, 2));
   
   const double beta = 1.5;
   
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
      
      // density
      db = (density_average[c](0) - density_average[c-1](0)) / dx;
      df = (density_average[c+1](0) - density_average[c](0)) / dx;
      dl = minmod ( density_average[c](1), beta * db, beta * df);
      
      density(local_dof_indices[0]) = density_average[c](0) - 
         0.5 * dx * dl;
      density(local_dof_indices[1]) = density_average[c](0) + 
         0.5 * dx * dl;
      
      // momentum
      db = (momentum_average[c](0) - momentum_average[c-1](0)) / dx;
      df = (momentum_average[c+1](0) - momentum_average[c](0)) / dx;
      dl = minmod ( momentum_average[c](1), beta * db, beta * df);
      
      momentum(local_dof_indices[0]) = momentum_average[c](0) - 
         0.5 * dx * dl;
      momentum(local_dof_indices[1]) = momentum_average[c](0) + 
         0.5 * dx * dl;
      
      // energy
      db = (energy_average[c](0) - energy_average[c-1](0)) / dx;
      df = (energy_average[c+1](0) - energy_average[c](0)) / dx;
      dl = minmod ( energy_average[c](1), beta * db, beta * df);
      
      energy(local_dof_indices[0]) = energy_average[c](0) - 
         0.5 * dx * dl;
      energy(local_dof_indices[1]) = energy_average[c](0) + 
         0.5 * dx * dl;
   }
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::compute_dt ()
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
   
   dt = std::min(dt, dx*dx/(mu_ref+1.0e-14));
   dt *= cfl;
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::output_results () const
{
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

    std::ofstream output ("solution.gpl");
    data_out.write_gnuplot (output);
   
   // compute shear stress and heat flux
   QIterated<dim>  quadrature = (fe.degree==1) ?
                                QIterated<dim>(QMidpoint<1>(),1)
                                :
                                QIterated<dim>(QTrapez<1>(),fe.degree-1);

   FEValues<dim> fe_values (fe, quadrature, update_quadrature_points |
                            update_values   | update_gradients);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points = quadrature.size();
   std::vector<double> density_values (n_q_points);
   std::vector<double> momentum_values (n_q_points);
   std::vector<double> energy_values (n_q_points);
   std::vector< Tensor<1,dim> >  density_grad (n_q_points);
   std::vector< Tensor<1,dim> >  momentum_grad (n_q_points);
   std::vector< Tensor<1,dim> >  energy_grad (n_q_points);

   
   typename DoFHandler<dim>::active_cell_iterator 
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();

   Vector<double> state (n_var);
   Vector<double> grad (n_var);
   double tau, q;
   
   std::ofstream fo;
   fo.open ("ns.dat");
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      fe_values.reinit (cell);
      fe_values.get_function_values (density, density_values);
      fe_values.get_function_values (momentum, momentum_values);
      fe_values.get_function_values (energy, energy_values);
      
      fe_values.get_function_gradients (density,  density_grad);
      fe_values.get_function_gradients (momentum, momentum_grad);
      fe_values.get_function_gradients (energy,   energy_grad);
      
      for(unsigned int i=0; i<n_q_points; ++i)
      {
         state(0) = density_values [i];
         state(1) = momentum_values[i];
         state(2) = energy_values  [i];
         
         grad(0) = density_grad [i][0];
         grad(1) = momentum_grad[i][0];
         grad(2) = energy_grad  [i][0];
         
         compute_nsterms (state, grad, tau, q);
         
         double pressure = (gas_gamma-1.0)*(energy_values[i] - 
                                            0.5 * pow(momentum_values[i],2) / density_values[i]);
         
         fo << fe_values.quadrature_point (i)(0) << " " 
            << tau << "  " << q << endl;
      }
      fo << endl;

   }
   
   fo.close ();
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::run ()
{
    std::cout << "Solving 1-D NS problem ...\n";

    make_grid_and_dofs();
    assemble_mass_matrix ();
    initialize ();
    output_results ();
    compute_averages ();

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
         //apply_limiter ();
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
       if(iter % 100 == 0) output_results ();
       
      std::cout << "Iter = " << iter << " time = " << time 
                << " Res =" << residual[0] << " " << residual[1] << " "
                << residual[2] << endl;
    }
    output_results ();
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main ()
{
    deallog.depth_console (0);
    {
        NSProblem<1> ns_problem(2, shock_structure);
        ns_problem.run ();
    }

    return 0;
}

