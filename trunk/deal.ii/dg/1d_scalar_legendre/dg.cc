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

#include <numerics/data_out.h>
#include <numerics/fe_field_function.h>

#include <base/convergence_table.h>

#include <fstream>
#include <iostream>

#include <base/logstream.h>

#define sign(a)   ((a>0.0) ? 1 : -1)

using namespace dealii;

// Number of variables: mass, momentum and energy
double xmin, xmax, xmid;
double u_left, u_right;
double final_time;
double Mlim; // used in TVB limiter
double Mh2;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[3] = {1.0, 1.0/4.0, 2.0/3.0};

// Numerical flux functions
enum FluxType {central, upwind};
enum TestCase {sine, hat};
enum LimiterType {none, tvd};

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
   LimiterType limiter_type;
   FluxType flux_type;
};

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a, const double& b, const double& c)
{
   double aa = std::fabs(a);
   if(aa < Mh2) return a;
   
   int sa = sign(a);
   int sb = sign(b);
   int sc = sign(c);
   
   double result;

   if( sa != sb || sb != sc)
   {
      result = 0.0;
   }
   else
   {
      result  = sa * std::min( aa, std::min(std::fabs(b), std::fabs(c)));
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
   InitialCondition (TestCase test_case)
   :
   Function<dim>(),
   test_case (test_case)
   {}
   
   virtual void value (const Point<dim>   &p,
                       double& values) const;
   
private:
   TestCase test_case;
};

// Initial condition
template<int dim>
void InitialCondition<dim>::value (const Point<dim>   &p,
                                   double& values) const
{
   double x = p[0];
   
   // test case: sine
   if(test_case == sine)
   {
      values = -std::sin (M_PI * x);
   }
   else if(test_case == hat)
   {
      if(std::fabs(x) < 0.25)
         values = 1.0;
      else
         values = 0.0;
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown test case"));
   }
}

//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class Solution : public Function<dim>
{
public:
   Solution (TestCase test_case)
   :
   Function<dim>(),
   test_case (test_case)
   {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
   virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                   const unsigned int  component = 0) const;
private:
   TestCase test_case;
};

// Exact solution works correctly only for periodic case
template<int dim>
double Solution<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
{
   double values;
   Point<dim> pp(p);
   pp[0] -= final_time;
   InitialCondition<dim> initial_condition(test_case);
   initial_condition.value(pp, values);
   return values;
}


// Exact solution works correctly only for periodic case
template<int dim>
Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
                                       const unsigned int) const
{
   Tensor<1,dim> values;
   double x = p[0] - final_time;
   
   // test case: sine
   if(test_case == sine)
   {
      values[0] = -M_PI * std::cos (M_PI * x);
   }
   else if(test_case == hat)
   {
      values[0] = 0;
   }
   else
   {
      AssertThrow(false, ExcMessage("Unknown test case"));
   }
   
   return values;
}



//------------------------------------------------------------------------------
// Flux of the PDE model
//------------------------------------------------------------------------------
double physical_flux (const double& u)
{
   return u;
}


//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
void CentralFlux (const double& left_state,
                  const double& right_state,
                  double& flux)
{
   flux = 0.5*(left_state + right_state);
}

//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
void UpwindFlux (const double& left_state,
                 const double& right_state,
                 double& flux)
{
   flux = left_state;
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
      case central:
         CentralFlux (left_state, right_state, flux);
         break;
         
      case upwind:
         UpwindFlux (left_state, right_state, flux);
         break;
         
      default:
         std::cout << "Unknown flux_type !!!\n";
         abort ();
   }
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
   void apply_TVD_limiter ();
   void mark_troubled_cells ();
   void update (const unsigned int rk_stage);
   void output_results (const double& time) const;
   void process_solution (unsigned int step);
   
   TestCase             test_case;
   bool                 debug;
   unsigned int         n_cells;
   double               dt;
   double               dx;
   double               cfl;
   unsigned int         n_rk_stages;
   LimiterType          limiter_type;
   FluxType             flux_type;
   unsigned int         nstep;
   
   
   Triangulation<dim>   triangulation;
   FE_DGP<dim>          fe;
   DoFHandler<dim>      dof_handler;
   
   std::vector< Vector<double> > inv_mass_matrix;
   
   Vector<double>       solution;
   Vector<double>       solution_old;
   Vector<double>       rhs;
   
   Vector<double>    average;
   std::vector<bool> troubled_cell;
   
   double residual;
   double residual0;
   
   ConvergenceTable  convergence_table;
   typename DoFHandler<dim>::active_cell_iterator firstc, lastc;
   std::vector<typename DoFHandler<dim>::active_cell_iterator> lcell, rcell;
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
    limiter_type (param.limiter_type),
    flux_type (param.flux_type),
    nstep (param.nstep),
    fe (param.degree),
    dof_handler (triangulation)
{
   AssertThrow (dim==1, ExcIndexRange(dim, 0, 1));
   
   final_time = param.final_time;
   
   n_rk_stages = 3;
   
   if(test_case == sine)
   {
      xmin    = -1.0;
      xmax    = +1.0;
   }
   else if(test_case == hat)
   {
      xmin    = -1.0;
      xmax    = +1.0;
   }
   else
   {
      std::cout << "Unknown test case\n";
   }
   
   dx = (xmax - xmin) / n_cells;
   Mh2 = Mlim*dx*dx;

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
   
    // allocate memory for inverse mass matrix. We store only diagonals.
    inv_mass_matrix.resize(n_cells);
    for (unsigned int c=0; c<n_cells; ++c)
       inv_mass_matrix[c].reinit(fe.degree+1);

   
   // Solution variables
   solution.reinit (dof_handler.n_dofs());
   solution_old.reinit (dof_handler.n_dofs());
   rhs.reinit (dof_handler.n_dofs());
   
   average.reinit (triangulation.n_cells());
   troubled_cell.resize(n_cells);
   
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
   
   Vector<double>   cell_matrix (dofs_per_cell);
   
   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
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
   
   InitialCondition<dim> initial_condition (test_case);
   
   double initial_value;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (unsigned int c=0; cell!=endc; ++cell, ++c)
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
            cell_rhs(i) += fe_values.shape_value (i, q_point) *
                           initial_value *
                           fe_values.JxW (q_point);
         }
      }
      
      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices (local_dof_indices);
      unsigned int ig, jg;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         ig = local_dof_indices[i];
         solution (ig) = cell_rhs(i) * inv_mass_matrix[c](i);
      }
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
            rhs (ig) = cell_rhs(i) * inv_mass_matrix[c](i);
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
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      cell->get_dof_indices (local_dof_indices);
      average[c] = solution(local_dof_indices[0]);
   }
}

//------------------------------------------------------------------------------
// Use TVB limiter to identify troubled cells
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::mark_troubled_cells ()
{
   const double EPS = 1.0e-14;
      
   QTrapez<dim>  quadrature_formula;
   
   FEValues<dim> fe_values (fe, quadrature_formula, update_values);
   std::vector<double> face_values(2), limited_face_values(2);
         
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   unsigned int n_troubled_cells = 0;
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(solution, face_values);

      unsigned int cl, cr;
      if(c==0)
      {
         cl = n_cells-1;
         cr = c+1;
      }
      else if(c==n_cells-1)
      {
         cl = c-1;
         cr = 0;
      }
      else
      {
         cl = c - 1;
         cr = c + 1;
      }
      double db = average[c]  - average[cl];
      double df = average[cr] - average[c];
      
      double DF = face_values[1] - average[c];
      double DB = average[c] - face_values[0];
      
      double dl = minmod ( DB, db, df);
      double dr = minmod ( DF, db, df);
      
      limited_face_values[0] = average[c] - dl;
      limited_face_values[1] = average[c] + dr;
      bool c0 = std::fabs(limited_face_values[0]-face_values[0])
                > EPS * std::fabs(face_values[0]);
      bool c1 = std::fabs(limited_face_values[1]-face_values[1])
                > EPS * std::fabs(face_values[1]);

      troubled_cell[c] = false;
      if(c0 || c1)
      {
         troubled_cell[c] = true;
         ++n_troubled_cells;
      }
   }
   //std::cout << "No of troubled cells = " << n_troubled_cells << std::endl;
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::apply_TVD_limiter ()
{
   if (fe.degree==0) return;
   
   static const double sqrt_3 = std::sqrt(3.0);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (unsigned int c=0; c<n_cells; ++c, ++cell)
   {
      cell->get_dof_indices (local_dof_indices);

      unsigned int cl, cr;
      if(c==0)
      {
         cl = n_cells-1;
         cr = c+1;
      }
      else if(c==n_cells-1)
      {
         cl = c-1;
         cr = 0;
      }
      else
      {
         cl = c - 1;
         cr = c + 1;
      }
      double db = average[c]  - average[cl];
      double df = average[cr] - average[c];
      
      double Dx = solution(local_dof_indices[1]);
      double Dx_new = minmod ( sqrt_3 * Dx, db, df) / sqrt_3;
      
      if(std::fabs(Dx - Dx_new) > 1.0e-6)
      {
         solution(local_dof_indices[1]) = Dx_new;
         for(unsigned int i=2; i<dofs_per_cell; ++i)
            solution(local_dof_indices[i]) = 0;
      }
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::apply_limiter ()
{
   if (fe.degree == 0 || limiter_type == none) return;
   apply_TVD_limiter ();
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
   std::cout << filename << std::endl;
   
   std::ofstream output (filename);
   data_out.write_gnuplot (output);
      
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   std::ofstream fo;
   filename = "avg_" + Utilities::int_to_string(c) + ".gpl";
   fo.open (filename);
   
   for (unsigned int c=0; cell!=endc; ++c, ++cell)
   {
      Point<dim> x = cell->center();
      int tc = (troubled_cell[c] ? 1 : 0);
      fo << x(0) << " " << average[c] << "  " << tc << std::endl;
   }
   
   fo.close ();
   
   ++c;
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
         apply_limiter ();
       }
       
       if(iter==0)
       {
          std::cout << "Initial residual = " << residual << std::endl;
          residual0 = residual;
       }
       
       residual /= residual0;
       
      time += dt;
      ++iter;
       if(iter % 100 == 0) output_results (time);
       
       if(debug)
      std::cout << "Iter = " << iter << " time = " << time 
                << " Res =" << residual << std::endl;
    }
   std::cout << "Iter = " << iter << " time = " << time 
             << " Res =" << residual << std::endl;

    output_results (time);
}
//------------------------------------------------------------------------------
// Compute error norms
//------------------------------------------------------------------------------
template <int dim>
void ScalarProblem<dim>::process_solution (unsigned int step)
{
   // compute error in solution
   Solution<dim> exact_solution(test_case);
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(2*fe.degree+1),
                                      VectorTools::L2_norm);
   const double L2_error = difference_per_cell.l2_norm();

   // compute error in gradient
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(2*fe.degree+1),
                                      VectorTools::H1_norm);
   const double H1_error = difference_per_cell.l2_norm();
   
   const unsigned int n_active_cells=triangulation.n_active_cells();
   const unsigned int n_dofs=dof_handler.n_dofs();
   convergence_table.add_value("step", step);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("L2", L2_error);
   convergence_table.add_value("H1", H1_error);
   
   convergence_table.set_scientific("L2", true);
   convergence_table.set_scientific("H1", true);
   
   std::cout << std::endl;
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
   convergence_table.evaluate_convergence_rates("L2",
                                                ConvergenceTable::reduction_rate_log2);

   convergence_table.set_precision("H1", 3);
   convergence_table.set_scientific("H1", true);
   convergence_table.evaluate_convergence_rates("H1",
                                                ConvergenceTable::reduction_rate_log2);

   convergence_table.set_tex_caption("cells", "\\# cells");
   convergence_table.set_tex_caption("dofs", "\\# dofs");
   convergence_table.set_tex_caption("L2", "$L^2$-error");
   convergence_table.set_tex_caption("H1", "$H^1$-error");
   
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
       param.degree       = 1;
       param.n_cells      = 50;
       param.nstep        = 1;
       param.test_case    = sine;
       param.cfl          = 0.9/(2.0*param.degree+1.0);
       param.final_time   = 10;
       param.limiter_type = none;
       param.flux_type    = upwind;
       
       bool debug = false;
       Mlim = 100.0;
       
       ScalarProblem<1> scalar_problem(param, debug);
       scalar_problem.run ();
    }

    return 0;
}

