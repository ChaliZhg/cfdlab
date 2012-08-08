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

// Number of variables: mass, sc and energy
const unsigned int n_var = 2;
double s_left, c_left, s_right, c_right;
double xmin, xmax, xmid;
unsigned int limit_what, limit_type;
double M_TVBM;
double beta;

const double velocity = 0.2;
const double gravity= 1.0;
const double density_oil = 1.0;
const double density_water = 2.0;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {dflu};
enum TestCase {t1};
//------------------------------------------------------------------------------

double viscosity_oil()
{
  return 1.0;
}
//------------------------------------------------------------------------------
double viscosity_water(double con)
{
   return 0.5 + con;
}
//------------------------------------------------------------------------------
double mobility_oil(double sat, double con)
{
   return std::pow(1.0-sat, 2.0) / viscosity_oil();
}
//------------------------------------------------------------------------------
double mobility_water(double sat, double con)
{
   return std::pow(sat,2.0) / viscosity_water(con);
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

// Initial condition for saturation and concentration
template<int dim>
void InitialCondition<dim>::vector_value (const Point<dim>   &p,
                                          Vector<double>& values) const
{

   if(p[0] < xmid)
   {
      values(0) = s_left; // left s
      values(1) = c_left; // left c
   }
   else
   {
      values(0) = s_right; // right s
      values(1) = c_right; // right c
   }

}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class ORProblem
{
public:
    ORProblem (unsigned int degree, unsigned int n_cells, TestCase test_case);
    void run ();

private:
    void make_grid_and_dofs ();
    void initialize ();
    void assemble_mass_matrix ();
    void assemble_rhs ();
    void compute_averages ();
    void compute_dt ();
    double minmod (const double& a, const double& b, const double& c);
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

    Triangulation<dim>   triangulation;
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> inv_mass_matrix;

    Vector<double>       s;
    Vector<double>       sc;
    Vector<double>       s_old;
    Vector<double>       sc_old;
    Vector<double>       rhs_s;
    Vector<double>       rhs_sc;
      
    std::vector< Vector<double> > s_average;
    std::vector< Vector<double> > sc_average;
   
    std::vector<double> residual;
    std::vector<double> residual0;

};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
ORProblem<dim>::ORProblem (unsigned int degree, unsigned int n_cells, TestCase test_case) :
    test_case (test_case),
    fe (degree),
    n_cells (n_cells),
    dof_handler (triangulation)
{
   Assert (dim==1, ExcIndexRange(dim, 0, 1));
      
   n_rk_stages = 3;
   flux_type = dflu;
   
   if(test_case == t1)
   {
      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.4;
      final_time = 1.0;
      min_residue= 1.0e20;
      cfl     = 0.4;
      
      s_left  = 0.1;
      s_right = 1.0;
      
      c_left  = 1.0;
      c_right = 0.0;
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
void ORProblem<dim>::make_grid_and_dofs ()
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
    s.reinit (dof_handler.n_dofs());
    s_old.reinit (dof_handler.n_dofs());
    rhs_s.reinit (dof_handler.n_dofs());
   
    sc.reinit (dof_handler.n_dofs());
    sc_old.reinit (dof_handler.n_dofs());
    rhs_sc.reinit (dof_handler.n_dofs());
   
    s_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
    sc_average.resize (triangulation.n_cells(), Vector<double>(fe.degree+1));
   
    residual.resize(n_var, 1.0);
    residual0.resize(n_var);
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::initialize ()
{
   std::cout << "Projecting initial condition ...";
   
   QGauss<dim>  quadrature_formula(fe.degree+1);
   
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   |
                            update_quadrature_points | 
                            update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   Vector<double>       cell_rhs_s  (dofs_per_cell);
   Vector<double>       cell_rhs_sc (dofs_per_cell);
   
   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   InitialCondition<dim> initial_condition;
   Vector<double> initial_value(n_var);
   double initial_s;
   double initial_sc;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      
      cell_rhs_s  = 0.0;
      cell_rhs_sc = 0.0;      
      
      // Flux integral over cell
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         // Get primitive variable at quadrature point
         initial_condition.vector_value(fe_values.quadrature_point(q_point),
                                        initial_value);
         // Convert primitive to conserved
         initial_s = initial_value(0);
         initial_sc= initial_value(0) * initial_value(1);

         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            cell_rhs_s(i) += (fe_values.shape_value (i, q_point) *
                              initial_s *
                              fe_values.JxW (q_point));
            cell_rhs_sc(i)+= (fe_values.shape_value (i, q_point) *
                              initial_sc *
                              fe_values.JxW (q_point));
         }
      }
      
      
      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices (local_dof_indices);
      unsigned int ig, jg;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         ig = local_dof_indices[i];
         
         s (ig) = 0.0;
         sc(ig) = 0.0;
         
         for (unsigned int j=0; j<dofs_per_cell; ++j)
         {
            jg = local_dof_indices[j];
            s(ig)  += inv_mass_matrix(ig,jg) * cell_rhs_s(j);
            sc(ig) += inv_mass_matrix(ig,jg) * cell_rhs_sc(j);
         }
         
      }
   }
   std::cout << "Done\n";
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::assemble_mass_matrix ()
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
// OR flux
//------------------------------------------------------------------------------
void or_flux (const double& sat,
              const double& sat_con,
              Vector<double>& flux)
{
  double con = (sat > 0) ? sat_con / sat : 0.0;
  double f=mobility_water(sat,con)/(mobility_water(sat,con)+mobility_oil(sat,con));
  flux(0)=(velocity-(density_water-density_oil)*gravity*mobility_oil(sat,con))*f;
  flux(1)=con*flux(0);
}
//----------------------------------------------------------------------------
// Find location of minimum for the flux
//----------------------------------------
double argmin_flux(const double& concentration)
{
   double r = viscosity_oil() / viscosity_water (concentration);
   double z = velocity * viscosity_oil() / 
      ((density_water - density_oil) * gravity);
   double alpha = 27.0 * (-r + r * r - z - 2.0 * r * z - r * r * z);
   double beta  = 2916.0 * r * r * r + alpha * alpha;

   double gamma = pow( alpha + sqrt(beta), 1.0/3.0);
   double fact  = 3.0 * pow(2.0, 1.0/3.0);

   double sstar = (1.0 - fact * r / gamma + gamma / fact) / (1.0 + r);

   double s_min;

   if(sstar <= 1.0 && sstar >= 0.0)
   {
      s_min = sstar;
   }
   else
      s_min = (velocity > 0.0) ? 0.0 : 1.0;

   return s_min;
}

//------------------------------------------------------------------------------
//----DFLU umerical flux function
//------------------------------------------------------------------------------
double DFLU_Flux
(
 const Vector<double>& left_state,
 const Vector<double>& right_state,
 Vector<double>& flux
 )
{
   double s_left  = left_state(0);
   double c_left  = left_state(1)/s_left;
   double s_right = right_state(0);
   double c_right = right_state(1)/s_right;
   
   if(gravity > 0.0) // Gravity is present
   {
      double s_min_left = argmin_flux(c_left);
      double s_min_right = argmin_flux(c_right);
      
      s_left  = std::max( s_left,  s_min_left);
      s_right = std::min( s_right, s_min_right);
      
      double m_water_left = mobility_water (s_left, c_left);
      double m_oil_left   = mobility_oil (s_left, c_left);
      double m_total_left = m_water_left + m_oil_left;
      
      double m_water_right = mobility_water (s_right, c_right);
      double m_oil_right   = mobility_oil (s_right, c_right);
      double m_total_right = m_water_right + m_oil_right;
      
      double v_left  = velocity 
      - (density_water - density_oil) * gravity * m_oil_left ;
      double v_right = velocity 
      - (density_water - density_oil) * gravity * m_oil_right ;
      double f_left  = v_left  * m_water_left / m_total_left;
      double f_right = v_right * m_water_right / m_total_right;
      flux(0)        = std::max( f_left, f_right );
      
      if(flux(0) > 0.0)
         flux(1) = c_left  * flux(0);
      else
         flux(1) = c_right * flux(0);
   }
   else // Gravity is not present
   {
      double m_water_left = mobility_water (s_left, c_left);
      double m_oil_left   = mobility_oil (s_left, c_left);
      double m_total_left = m_water_left + m_oil_left;
      
      double m_water_right = mobility_water (s_right, c_right);
      double m_oil_right   = mobility_oil (s_right, c_right);
      double m_total_right = m_water_right + m_oil_right;
      
      if (velocity > 0)
      {
         flux(0) = velocity * m_water_left / m_total_left;
         flux(1) = c_left * flux(0);
      }
      else
      {
         flux(0) = velocity * m_water_right / m_total_right;
         flux(1) = c_right * flux(0);
      }
   }
   
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
      case dflu:
         DFLU_Flux (left_state, right_state, flux);
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
void ORProblem<dim>::assemble_rhs ()
{
    QGaussLobatto<dim>  quadrature_formula(fe.degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients |
                             update_quadrature_points | 
                             update_JxW_values);

   // for getting neighbour cell solutions to compute intercell flux
   QTrapez<dim> quadrature_dummy;
   FEValues<dim> fe_values_neighbor (fe, quadrature_dummy,
                            update_values   | update_gradients);
   
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<double>  s_values  (n_q_points);
    std::vector<double>  sc_values (n_q_points);
   
   // for getting neighbor cell solution using trapezoidal rule
   std::vector<double>  s_values_n  (2);
   std::vector<double>  sc_values_n (2);
   
    Vector<double>       cell_rhs_s  (dofs_per_cell);
    Vector<double>       cell_rhs_sc (dofs_per_cell);
   
    Vector<double>       flux(n_var);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
    residual[0] = residual[1] = 0.0;
   
    for (unsigned int c=0; cell!=endc; ++cell, ++c)
    {
        fe_values.reinit (cell);
       
        cell_rhs_s  = 0.0;
        cell_rhs_sc = 0.0;

        // Compute conserved variables at quadrature points
        fe_values.get_function_values (s,  s_values);
        fe_values.get_function_values (sc, sc_values);

        // Flux integral over cell
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            or_flux(s_values[q_point], sc_values[q_point], flux);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                cell_rhs_s(i) += (fe_values.shape_grad (i, q_point)[0] *
                                  flux(0) *
                                  fe_values.JxW (q_point));
                cell_rhs_sc(i)+= (fe_values.shape_grad (i, q_point)[0] *
                                  flux(1) *
                                  fe_values.JxW (q_point));
            }
        }
       
       // Computation of flux at cell boundaries
       Vector<double> lf_left_state(n_var), lf_right_state(n_var);
       Vector<double> lf_left_grad(n_var), lf_right_grad(n_var);

        // left face flux
        // right state is from current cell
       lf_right_state(0) = s_values [0];
       lf_right_state(1) = sc_values[0];
       
       if(c==0)
       {
          lf_left_state(0) = s_left;
          lf_left_state(1) = s_left * c_left;
       }
       else
       {
          // get left cell dof indices
          fe_values_neighbor.reinit (cell->neighbor(0));
          
          fe_values_neighbor.get_function_values (s,  s_values_n);
          fe_values_neighbor.get_function_values (sc, sc_values_n);
          
          lf_left_state(0) = s_values_n [1];
          lf_left_state(1) = sc_values_n[1];
       }
       
       Vector<double> left_flux(n_var);
       numerical_flux (flux_type, lf_left_state, lf_right_state, left_flux);
       
       // right face flux
       Vector<double> rf_left_state(n_var), rf_right_state(n_var);
       // left state is from current cell
       rf_left_state(0) = s_values [n_q_points-1];
       rf_left_state(1) = sc_values[n_q_points-1];
       
       if(c==triangulation.n_cells()-1)
       {
          rf_right_state(0) = s_right;
          rf_right_state(1) = s_right * c_right;
       }
       else
       {          
          // get right cell to right face
          fe_values_neighbor.reinit (cell->neighbor(1));
          
          fe_values_neighbor.get_function_values (s,  s_values_n);
          fe_values_neighbor.get_function_values (sc, sc_values_n);
          
          rf_right_state(0) = s_values_n [0];
          rf_right_state(1) = sc_values_n[0];
       }
              
       Vector<double> right_flux(3);
       numerical_flux (flux_type, rf_left_state, rf_right_state, right_flux);
              
        // Add flux at cell boundaries
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
           // Left face flux
           cell_rhs_s(i) += fe_values.shape_value (i, 0) *
                            left_flux(0);
           cell_rhs_sc(i)+= fe_values.shape_value (i, 0) *
                            left_flux(1);
           
           // Right face flux
           cell_rhs_s(i) -= fe_values.shape_value (i, n_q_points-1) *
                            right_flux(0);
           cell_rhs_sc(i)-= fe_values.shape_value (i, n_q_points-1) *
                            right_flux(1);
        }

        // Multiply by inverse mass matrix and add to rhs
        cell->get_dof_indices (local_dof_indices);
        unsigned int ig, jg;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            ig = local_dof_indices[i];
           
            rhs_s (ig) = 0.0;
            rhs_sc(ig) = 0.0;
           
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
               jg = local_dof_indices[j];
               rhs_s(ig)  += inv_mass_matrix(ig,jg) * cell_rhs_s(j);
               rhs_sc(ig) += inv_mass_matrix(ig,jg) * cell_rhs_sc(j);
            }
           
            residual[0] += std::pow (rhs_s (ig), 2);
            residual[1] += std::pow (rhs_sc(ig), 2);
        }
       
    }

}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::compute_averages ()
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
      
      s_average[c]  = 0.0;
      sc_average[c] = 0.0;
      for(unsigned int point=0; point<n_q_points; ++point)
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            s_average[c](0) += s(local_dof_indices[i]) * 
                               fe_values.shape_value (i, point) *
                               fe_values.JxW (point);
            sc_average[c](0) += sc(local_dof_indices[i]) * 
                                fe_values.shape_value (i, point) *
                                fe_values.JxW (point);
            
            if(fe.degree >= 1) // compute average gradient
            {
               s_average[c](1) += s(local_dof_indices[i]) * 
                                  fe_values.shape_grad (i, point)[0] *
                                  fe_values.JxW (point);
               sc_average[c](1) += sc(local_dof_indices[i]) * 
                                   fe_values.shape_grad (i, point)[0] *
                                   fe_values.JxW (point);
            }
         }
      
      s_average[c]  /= dx;
      sc_average[c] /= dx;
   }
}

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
template <int dim>
double ORProblem<dim>::minmod (const double& a, const double& b, const double& c)
{
   double result;
   if(std::fabs(a) < M_TVBM * dx * dx)
   {
      result = a;
   }
   else if( a*b > 0.0 && b*c > 0.0)
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
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::apply_limiter ()
{
   Assert (fe.degree==1, ExcIndexRange(fe.degree, 1, 2));
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;   
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   // dont limit in first cell, skip it.
   ++cell;
   
   for (unsigned int c=1; c<n_cells-1; ++c, ++cell)
   {
      cell->get_dof_indices (local_dof_indices);
      
      // limit s
      double dbs = (s_average[c](0) - s_average[c-1](0)) / dx;
      double dfs = (s_average[c+1](0) - s_average[c](0)) / dx;
      double dls = minmod ( s_average[c](1), beta*dbs, beta*dfs);
      
      s(local_dof_indices[0]) = s_average[c](0) - 0.5 * dx * dls;
      s(local_dof_indices[1]) = s_average[c](0) + 0.5 * dx * dls;
      
      // limit sc
      if(limit_what == 0)
      {
         double dbsc = (sc_average[c](0) - sc_average[c-1](0)) / dx;
         double dfsc = (sc_average[c+1](0) - sc_average[c](0)) / dx;
         double dlsc = minmod ( sc_average[c](1), beta*dbsc, beta*dfsc);
         
         sc(local_dof_indices[0]) = sc_average[c](0) - 0.5 * dx * dlsc;
         sc(local_dof_indices[1]) = sc_average[c](0) + 0.5 * dx * dlsc;
      }
      else
      // Limit c and then get sc
      {
         double sl = s_average[c](0) - 0.5 * dx * s_average[c](1);
         double sr = s_average[c](0) + 0.5 * dx * s_average[c](1);
         double scl = sc_average[c](0) - 0.5 * dx * sc_average[c](1);
         double scr = sc_average[c](0) + 0.5 * dx * sc_average[c](1);

         double dbc = (sc_average[c](0)/s_average[c](0) - sc_average[c-1](0)/s_average[c-1](0)) / dx;
         double dfc = (sc_average[c+1](0)/s_average[c+1](0) - sc_average[c](0)/s_average[c](0)) / dx;
         double dcc = (scr/sr - scl/sl) / dx;
         double dlc = minmod ( dcc, beta*dbc, beta*dfc);
         
         // limited slope of sc
         double dlsc = sc_average[c](0)/s_average[c](0) * dls + s_average[c](0) * dlc;
         sc(local_dof_indices[0]) = sc_average[c](0) - 0.5 * dx * dlsc;
         sc(local_dof_indices[1]) = sc_average[c](0) + 0.5 * dx * dlsc;
      }
      
   }
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cells; ++i)
   {
      dt = dx/0.8;
   }
   
   dt *= cfl;
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::update (const unsigned int rk_stage)
{
   // Update conserved variables
   for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   {
      s(i)  = a_rk[rk_stage] * s_old(i) +
                    b_rk[rk_stage] * (s(i) + dt * rhs_s(i));
      sc(i) = a_rk[rk_stage] * sc_old(i) +
                    b_rk[rk_stage] * (sc(i) + dt * rhs_sc(i));
   }

}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::output_results () const
{
    Vector<double> c(dof_handler.n_dofs());

    // Compute velocity and pressure
    for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
       c(i) = sc(i) / s(i);
    }

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (s, "saturation");
    data_out.add_data_vector (c, "concentration");

    if(fe.degree <= 1)
       data_out.build_patches (1);
    else
       data_out.build_patches (fe.degree);

    std::ofstream output ("solution.gpl");
    data_out.write_gnuplot (output);
   
   // save cell averages
   std::ofstream avg("averages.gpl");
   for(unsigned int c=0; c<n_cells; ++c)
   {
      double con = sc_average[c](0) / s_average[c](0);
      double x = xmin + c*dx + 0.5*dx;
      avg << x << "  " << s_average[c](0) << "  " << con << std::endl;
   }
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void ORProblem<dim>::run ()
{
    std::cout << "Solving 1-D oil reservoir problem ...\n";

    make_grid_and_dofs();
    assemble_mass_matrix ();
    initialize ();
    output_results ();
    compute_averages ();
    apply_limiter ();

    double time = 0.0;
    unsigned int iter = 0;

    while (time < final_time)
    {
       s_old  = s;
       sc_old = sc;
       
       compute_dt ();
       if(time+dt > final_time) dt = final_time - time;

       for(unsigned int rk=0; rk<n_rk_stages; ++rk)
       {
         assemble_rhs ();
         update (rk);
         compute_averages ();
         apply_limiter ();
       }
       
       if(iter==0)
       {
          std::cout << "Initial residual = " << residual[0] << " "
                    << residual[1] << endl;
          for(unsigned int i=0; i<n_var; ++i)
             residual0[i] = residual[i];
       }
       
       for(unsigned int i=0; i<n_var; ++i)
          residual[i] /= residual0[i];
       
      time += dt;
      ++iter;
       if(iter % 100 == 0) output_results ();
       
      std::cout << "Iter = " << iter << " time = " << time 
                << " Res =" << residual[0] << " "
                << residual[1] << endl;
    }
    output_results ();
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main ()
{
   int degree;
   int n_cells;
   std::cout << "Degree of basis function = ";
   std::cin >> degree;
   std::cout << "Number of cells = ";
   std::cin >> n_cells;
   std::cout << "Limit (0=sc, 1=c) = ";
   std::cin >> limit_what;
   assert(limit_what==0 || limit_what==1);
   std::cout << "Limiter (0=TVDM, 1=TVBM) = ";
   std::cin >> limit_type;
   assert(limit_type==0 || limit_type==1);
   if(limit_type == 0)
   {
      std::cout << "   beta = ";
      std::cin >> beta;
      assert(beta > 0.0);
      M_TVBM = 0.0;
   }
   else if(limit_type == 1)
   {
      beta = 1.0;
      std::cout << "   M = ";
      std::cin >> M_TVBM;
      assert(M_TVBM > 0.0);
   }
   
   deallog.depth_console (0);
   {
      ORProblem<1> or_problem(degree, n_cells, t1);
      or_problem.run ();
   }
   
   return 0;
}
