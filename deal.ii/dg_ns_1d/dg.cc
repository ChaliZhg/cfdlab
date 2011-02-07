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

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {lxf, kfvs};
enum TestCase {sod, shock_structure};

//------------------------------------------------------------------------------
// Viscosity coefficient as function of temperature
//------------------------------------------------------------------------------
double viscosity (const double& T)
{
   return mu_ref * pow( T/T_ref, 0.8);
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
   else {
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
    void compute_face_flux ();
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
    double                cfl;
    double               final_time;
    unsigned int         n_rk_stages;
    FluxType             flux_type;


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
   

   n_cells = 100;

   n_rk_stages = 3;
   flux_type = kfvs;
   
   if(test_case == sod)
   {
      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.5;
      final_time = 2.0;
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
   }
   else if (test_case == shock_structure)
   {
      xmin    = -0.25;
      xmax    =  0.0;
      xmid    = 0.5 * ( xmin + xmax );
      final_time = 2.0;
      cfl = 0.1;
      
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
   }
   else
   {
      std::cout << "Unknown test case\n";
   }
   
   dx      = (xmax - xmin) / n_cells;

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
void kfvs_split_flux (const std::vector<double>& prim,
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
   
   // inviscid flux
   flux[0] = prim[0] * fact;
   flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A +
   prim[0] * prim[1] * B;
   flux[2] = prim[1] * (E + prim[2]) * A +
             (E + 0.5 * prim[2]) * B;
   
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
// KFVS flux
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
   
   kfvs_split_flux (left,  tau_l, q_l, +1, flux_pos);
   kfvs_split_flux (right, tau_r, q_r, -1, flux_neg);
   
   for(unsigned int i=0; i<n_var; ++i)
      flux(i) = flux_pos[i] + flux_neg[i];
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
template <int dim>
void NSProblem<dim>::compute_face_flux ()
{
   const unsigned int dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> dofs(dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      l_cell, 
      r_cell;
   
   unsigned int l_dof, r_dof;
   
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
         
      Vector<double> left_state(n_var), right_state(n_var);

      left_state(0) = density  (l_dof);
      left_state(1) = momentum (l_dof);
      left_state(2) = energy   (l_dof);
      
      right_state(0) = density  (r_dof);
      right_state(1) = momentum (r_dof);
      right_state(2) = energy   (r_dof);
      
      
      std::vector< Tensor<1,dim> > l_grad(n_var);
      std::vector< Tensor<1,dim> > r_grad(n_var);
      
      for(unsigned int k=0; k<n_var; ++k)
      {
         l_grad[k][0] = 0.0;
         r_grad[k][0] = 0.0;
      }
            
      // This part of code is VERY SLOW
      if(test_case == shock_structure)
      {
         Functions::FEFieldFunction<dim> density_grad (dof_handler, density);
         density_grad.set_active_cell (l_cell);
         l_grad[0] = density_grad.gradient (l_cell->vertex(1));
         density_grad.set_active_cell (r_cell);
         r_grad[0] = density_grad.gradient (r_cell->vertex(0));

         Functions::FEFieldFunction<dim> momentum_grad (dof_handler, momentum);
         momentum_grad.set_active_cell (l_cell);
         l_grad[1] = momentum_grad.gradient (l_cell->vertex(1));
         momentum_grad.set_active_cell (r_cell);
         r_grad[1] = momentum_grad.gradient (r_cell->vertex(0));

         Functions::FEFieldFunction<dim> energy_grad (dof_handler, energy);
         energy_grad.set_active_cell (l_cell);
         l_grad[2] = energy_grad.gradient (l_cell->vertex(1));
         energy_grad.set_active_cell (r_cell);
         r_grad[2] = energy_grad.gradient (r_cell->vertex(0));
      }
      
      Vector<double> left_grad (n_var), right_grad (n_var);
      for(unsigned int k=0; k<n_var; ++k)
      {
         left_grad(k) = l_grad[k][0];
         right_grad(k) = r_grad[k][0];
      }

      
      switch (flux_type) {
         case lxf:
            LaxFlux (left_state, right_state, face_flux[i]);
            break;
            
         case kfvs:
            KFVSFlux (left_state, right_state, left_grad, right_grad, face_flux[i]);
            break;
            
         default:
            cout << "Unknown flux_type !!!\n";
            abort ();
      }
   }
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

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  density_values  (n_q_points);
    std::vector<double>  momentum_values (n_q_points);
    std::vector<double>  energy_values   (n_q_points);
   
    std::vector< Tensor<1,dim> >  density_grad (n_q_points);
    std::vector< Tensor<1,dim> >  momentum_grad (n_q_points);
    std::vector< Tensor<1,dim> >  energy_grad (n_q_points);


    Vector<double>       cell_rhs_density  (dofs_per_cell);
    Vector<double>       cell_rhs_momentum (dofs_per_cell);
    Vector<double>       cell_rhs_energy   (dofs_per_cell);
   
    Vector<double>       flux(n_var);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
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

        // Add flux at cell boundaries
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

    data_out.build_patches ();

    std::ofstream output ("solution.gpl");
    data_out.write_gnuplot (output);
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

    while (time < final_time)
    {
       density_old  = density;
       momentum_old = momentum;
       energy_old   = energy;
       
       compute_dt ();
       if(time+dt > final_time) dt = final_time - time;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
         compute_face_flux ();
         assemble_rhs ();
         update (rk);
         compute_averages ();
         apply_limiter ();
       }
      time += dt;
      ++iter;
      std::cout << "Iter = " << iter << " time = " << time << endl;
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
        NSProblem<1> ns_problem(1, shock_structure);
        ns_problem.run ();
    }

    return 0;
}

