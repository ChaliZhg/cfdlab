/* $Id: step-12.cc 22886 2010-11-29 23:48:38Z bangerth $ */
/* Author: Guido Kanschat, Texas A&M University, 2009 */

/*    $Id: step-12.cc 22886 2010-11-29 23:48:38Z bangerth $       */
/*                                                                */
/*    Copyright (C) 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

// Modifications by Praveen. C, http://math.tifrbng.res.in/~praveen
// Explicit time-stepping Runge-Kutta DG method
// Mass matrix on each cell is computed, inverted and the inverse 
// is stored. Then in each time iteration, we need to compute right
// hand side and multipy by inverse mass mass matrix. After that
// solution is advanced to new time level by an RK scheme.

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <fe/fe_values.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/data_out.h>
#include <numerics/vector_tools.h>
#include <fe/mapping_q1.h>
		
#include <fe/fe_dgp.h>
				 
#include <numerics/derivative_approximation.h>
#include <numerics/solution_transfer.h>

#include <meshworker/dof_info.h>
#include <meshworker/integration_info.h>
#include <meshworker/simple.h>
#include <meshworker/loop.h>

#include <iostream>
#include <fstream>
#include <cmath>

#define sign(a)  ( (a) > 0 ? 1 : -1 )

using namespace dealii;

const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};
double Mlim;
double dx;

enum LimiterType {none, tvd};
enum TestCase {expo, square, circ};

//------------------------------------------------------------------------------
// Minmod function
//------------------------------------------------------------------------------
double minmod(double a, double b, double c)
{

   double aa = std::fabs(a);
   if(aa < Mlim*dx*dx) return a;

   int sa = sign(a);
   int sb = sign(b);
   int sc = sign(c);

   if(sa != sb) return 0;
   if(sa != sc) return 0;

   return sa * std::min( aa, std::min( sb*b, sc*c ) );
}

//------------------------------------------------------------------------------
// Speed of advection
//------------------------------------------------------------------------------
template <int dim>
void advection_speed(const Point<dim>& p, Point<dim>& v)
{
   v(0) = -p(1);
   v(1) =  p(0);
}
//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
double numerical_flux(double vel, double ul, double ur)
{
   if(vel > 0)
      return vel * ul;
   else
      return vel * ur;
}
//------------------------------------------------------------------------------
// Initial condition function class
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition: public Function<dim>
{
public:
   InitialCondition () {};
   virtual void value_list (const std::vector<Point<dim> > &points,
                            std::vector<double> &values,
                            const unsigned int component=0) const;
   TestCase test_case;
};

// Computes boundary condition value at a list of boundary points
template <int dim>
void InitialCondition<dim>::value_list(const std::vector<Point<dim> > &points,
                                     std::vector<double> &values,
                                     const unsigned int) const
{
   Assert(values.size()==points.size(),
          ExcDimensionMismatch(values.size(),points.size()));
   
   if(test_case == expo)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         double r2 = std::pow(points[i](0)-0.5, 2.0) + std::pow(points[i](1), 2.0);
         double r = std::sqrt(r2);
         if(r < 0.25)
            values[i] = 1.0 + std::pow(std::cos(2.0*M_PI*r), 2.0);
         else
            values[i] = 1.0;
      }
   else if(test_case == circ)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         double r2 = std::pow(points[i](0)-0.5, 2.0) + std::pow(points[i](1), 2.0);
         if(r2 < 0.25*0.25)
            values[i] = 2.0;
         else
            values[i] = 1.0;
      }
   else if(test_case == square)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         const Point<dim>& p = points[i];
         if(std::fabs(p(0)-0.5) < 0.25 && std::fabs(p(1)) < 0.25)
            values[i] = 2.0;
         else
            values[i] = 1.0;
      }
   else
      for (unsigned int i=0; i<values.size(); ++i)
         values[i] = 1.0;
}

//------------------------------------------------------------------------------
// Boundary condition function class
//------------------------------------------------------------------------------
template <int dim>
class BoundaryValues: public Function<dim>
{
  public:
    BoundaryValues () {};
    virtual void value_list (const std::vector<Point<dim> > &points,
			                    std::vector<double> &values,
			                    const unsigned int component=0) const;
};

// Computes boundary condition value at a list of boundary points
template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim> > &points,
				       std::vector<double> &values,
				       const unsigned int) const
{
   Assert(values.size()==points.size(),
          ExcDimensionMismatch(values.size(),points.size()));
   
   for (unsigned int i=0; i<values.size(); ++i)
   {
      values[i]=1.0;
   }
}

//------------------------------------------------------------------------------
// Class for integrating rhs using MeshWorker
//------------------------------------------------------------------------------
template <int dim>
class RHSIntegrator
{
   public:
      RHSIntegrator (const DoFHandler<dim>& dof_handler)
         : dof_info (dof_handler) {};

      MeshWorker::IntegrationInfoBox<dim> info_box;
      MeshWorker::DoFInfo<dim> dof_info;
      MeshWorker::Assembler::ResidualSimple< Vector<double> >
         assembler;
};

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class Step12
{
   public:
      Step12 (unsigned int degree, 
              LimiterType  limiter_type,
              TestCase     test_case);
      void run ();
      
   private:
      void setup_system ();
      void assemble_mass_matrix ();
      void set_initial_condition ();
      void setup_mesh_worker (RHSIntegrator<dim>&);
      void assemble_rhs (RHSIntegrator<dim>&);
      void compute_dt ();
      void solve ();
      void compute_shock_indicator ();
      void apply_limiter_old ();
      void apply_limiter ();
      void apply_limiter_TVD ();
      void refine_grid ();
      void compute_min_max ();
      void output_results (const unsigned int cycle);
   double compute_cell_average(const typename DoFHandler<dim>::cell_iterator& cell);
      
      Triangulation<dim>   triangulation;
      const MappingQ1<dim> mapping;
   
      unsigned int         degree;
      FE_DGP<dim>          fe;
      DoFHandler<dim>      dof_handler;
      FE_DGP<dim>          fe_cell;
      DoFHandler<dim>      dof_handler_cell;
   
      std::vector< Vector<double> > inv_mass_matrix;
      
      Vector<double>       solution;
      Vector<double>       solution_old;
      Vector<double>       average;
      Vector<double>       right_hand_side;
      Vector<double>       shock_indicator;
      Vector<double>       jump_indicator;
      double               dt;
      double               cfl;
      LimiterType          limiter_type;
      TestCase             test_case;
      double               sol_min, sol_max;
   
      std::vector<typename DoFHandler<dim>::cell_iterator>
         lcell, rcell, bcell, tcell;
   
      typedef MeshWorker::DoFInfo<dim> DoFInfo;
      typedef MeshWorker::IntegrationInfo<dim> CellInfo;
      
      static void integrate_cell_term (DoFInfo& dinfo, CellInfo& info);
      static void integrate_boundary_term (DoFInfo& dinfo, CellInfo& info);
      static void integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                       CellInfo& info1, CellInfo& info2);
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
Step12<dim>::Step12 (unsigned int degree,
                     LimiterType  limiter_type,
                     TestCase     test_case)
      :
      mapping (),
      degree (degree),
      fe (degree),
      dof_handler (triangulation),
      fe_cell(0),
      dof_handler_cell (triangulation),
      limiter_type (limiter_type),
      test_case (test_case)
{
   cfl = 0.8/(2.0*degree + 1.0);
}

//------------------------------------------------------------------------------
// Make dofs and allocate memory
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_system ()
{
   std::cout << "Allocating memory ...\n";

   dof_handler.distribute_dofs (fe);
   dof_handler_cell.distribute_dofs (fe_cell);

   inv_mass_matrix.resize(triangulation.n_cells());
   for (unsigned int c=0; c<triangulation.n_cells(); ++c)
      inv_mass_matrix[c].reinit(fe.dofs_per_cell);

   solution.reinit (dof_handler.n_dofs());
   solution_old.reinit (dof_handler.n_dofs());
   right_hand_side.reinit (dof_handler.n_dofs());
   
   average.reinit (dof_handler_cell.n_dofs());
   shock_indicator.reinit (dof_handler_cell.n_dofs());
   jump_indicator.reinit (dof_handler_cell.n_dofs());

   // For each cell, find neighbourig cell
   // This is needed for limiter
   lcell.resize(triangulation.n_cells());
   rcell.resize(triangulation.n_cells());
   bcell.resize(triangulation.n_cells());
   tcell.resize(triangulation.n_cells());
   
   const double EPS = 1.0e-10;
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      lcell[c] = endc;
      rcell[c] = endc;
      bcell[c] = endc;
      tcell[c] = endc;
      dx = cell->diameter() / std::sqrt(2.0);
      
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
         if (! cell->at_boundary(face_no))
         {
            const typename DoFHandler<dim>::cell_iterator
            neighbor = cell->neighbor(face_no);
            Assert(neighbor->level() == cell->level() || neighbor->level() == cell->level()-1,
                   ExcInternalError());
            Point<dim> dr = neighbor->center() - cell->center();
            if(dr(0) < -0.5*dx)
               lcell[c] = neighbor;
            else if(dr(0) > 0.5*dx)
               rcell[c] = neighbor;
            else if(dr(1) < -0.5*dx)
               bcell[c] = neighbor;
            else if(dr(1) > 0.5*dx)
               tcell[c] = neighbor;
            else
            {
               std::cout << "Did not find all neighbours\n";
               std::cout << "dx, dy = " << dr(0) << "  " << dr(1) << std::endl;
               exit(0);
            }
         }
   }
   
   assemble_mass_matrix ();
   
   std::cout << "Number of active cells:       "
             << triangulation.n_active_cells()
             << std::endl;
   
   std::cout << "Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << std::endl;
}
//------------------------------------------------------------------------------
// If cell is active, return cell average.
// If is not active, return arithmetic average of child cells.
//------------------------------------------------------------------------------
template <int dim>
double Step12<dim>::compute_cell_average(const typename DoFHandler<dim>::cell_iterator& cell)
{
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   
   if(cell->active())
   {
      cell->get_dof_indices(dof_indices);
      return solution(dof_indices[0]);
   }
   else
   {  // compute average solution on child cells
      auto child_cells = GridTools::get_active_child_cells< DoFHandler<dim> > (cell);
      double avg = 0, measure = 0;
      for(unsigned int i=0; i<child_cells.size(); ++i)
      {
         child_cells[i]->get_dof_indices(dof_indices);
         avg += solution(dof_indices[0]) * child_cells[i]->measure();
         measure += child_cells[i]->measure();
      }
      avg /= measure;
      return avg;
   }
}
//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_mass_matrix ()
{
   std::cout << "Constructing mass matrix ...\n";
   
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
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      fe_values.reinit (cell);
      cell_matrix = 0.0;
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
               cell_matrix(i) += fe_values.shape_value (i, q_point) *
                                 fe_values.shape_value (i, q_point) *
                                 fe_values.JxW (q_point);
      
      // Invert cell_matrix
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         inv_mass_matrix[c](i) = 1.0/cell_matrix(i);
   }
   
}
//------------------------------------------------------------------------------
// Project initial condition
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::set_initial_condition ()
{
   InitialCondition<dim> initial_condition;
   initial_condition.test_case = test_case;
   VectorTools::create_right_hand_side(dof_handler,
                                       QGauss<dim>(fe.degree+1),
                                       initial_condition,
                                       solution);
   
   // Multiply by inverse mass matrix
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      cell->get_dof_indices (local_dof_indices);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         solution(local_dof_indices[i]) *= inv_mass_matrix[c](i);
   }

   apply_limiter ();
}
//------------------------------------------------------------------------------
// Create mesh worker for integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_mesh_worker (RHSIntegrator<dim>& rhs_integrator)
{   
   std::cout << "Setting up mesh worker ...\n";

   MeshWorker::IntegrationInfoBox<dim>& info_box = rhs_integrator.info_box;
   MeshWorker::DoFInfo<dim>& dof_info = rhs_integrator.dof_info;
   MeshWorker::Assembler::ResidualSimple< Vector<double> >&
      assembler = rhs_integrator.assembler;

   const unsigned int n_gauss_points = fe.degree+1;
   info_box.initialize_gauss_quadrature(n_gauss_points,
                                        n_gauss_points,
                                        n_gauss_points);

   // Add solution vector to info_box
   NamedData< Vector<double>* > solution_data;
   solution_data.add (&solution, "solution");
   info_box.cell_selector.add     ("solution", true, false, false);
   info_box.boundary_selector.add ("solution", true, false, false);
   info_box.face_selector.add     ("solution", true, false, false);
   
   info_box.initialize_update_flags ();
   info_box.add_update_flags_all      (update_quadrature_points);
   info_box.add_update_flags_cell     (update_gradients);
   info_box.add_update_flags_boundary (update_values);
   info_box.add_update_flags_face     (update_values);

   info_box.initialize (fe, mapping, solution_data);
   
   // Attach rhs vector to assembler
   NamedData< Vector<double>* > rhs;
   Vector<double>* data = &right_hand_side;
   rhs.add (data, "RHS");
   assembler.initialize (rhs);
}

//------------------------------------------------------------------------------
// Compute time-step
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::compute_dt ()
{
   std::cout << "Computing local time-step ...\n";
      
   dt = 1.0e20;
   
   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      double h = cell->diameter ();
      const Point<dim> cell_center = cell->center();
      Point<dim> beta;
      advection_speed(cell_center, beta);
      
      dt = std::min ( dt, h / beta.norm ());
   }
   
   dt *= cfl;
   
}

//------------------------------------------------------------------------------
// Assemble rhs of the problem
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_rhs (RHSIntegrator<dim>& rhs_integrator)
{
   right_hand_side = 0.0;

   MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
                    MeshWorker::IntegrationInfoBox<dim> >
      (dof_handler.begin_active(), dof_handler.end(),
       rhs_integrator.dof_info, 
       rhs_integrator.info_box,
       &Step12<dim>::integrate_cell_term,
       &Step12<dim>::integrate_boundary_term,
       &Step12<dim>::integrate_face_term,
       rhs_integrator.assembler, true);

   // Multiply by inverse mass matrix
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   {
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
         right_hand_side(local_dof_indices[i]) *= inv_mass_matrix[c](i);
   }
}

//------------------------------------------------------------------------------
// Compute cell integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_cell_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);
      
      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
         local_vector(i) += beta * fe_v.shape_grad(i,point) *
                            sol[point] * JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute boundary integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_boundary_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];

   Vector<double>& local_vector = dinfo.vector(0).block(0);
   
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Point<dim> >& normals = fe_v.get_normal_vectors ();
   
   std::vector<double> g(fe_v.n_quadrature_points);
   
   static BoundaryValues<dim> boundary_function;
   boundary_function.value_list (fe_v.get_quadrature_points(), g);
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);
      
      const double beta_n = beta * normals[point];
      const double flux = numerical_flux(beta_n, sol[point], g[point]);
      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
         local_vector(i) -= flux *
                            fe_v.shape_value(i,point) *
                            JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute integral over internal faces
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
				                           CellInfo& info1, CellInfo& info2)
{
   const FEValuesBase<dim>& fe_v          = info1.fe_values();
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();
   
   const std::vector<double>& sol1 = info1.values[0][0];
   const std::vector<double>& sol2 = info2.values[0][0];
   
   Vector<double>& local_vector1 = dinfo1.vector(0).block(0);
   Vector<double>& local_vector2 = dinfo2.vector(0).block(0);
   
   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Point<dim> >& normals = fe_v.get_normal_vectors ();
   
   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);
      
      const double beta_n = beta * normals[point];
      const double flux = numerical_flux(beta_n, sol1[point], sol2[point]);
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= flux *
                                fe_v.shape_value(i,point) *
                                JxW[point];
         
         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += flux *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
   }
}
//------------------------------------------------------------------------------
// Compute KXRCF indicator
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::compute_shock_indicator ()
{
   QGauss<dim-1> quadrature(fe.degree + 1);
   FEFaceValues<dim> fe_face_values (fe, quadrature,
                                     update_values | update_normal_vectors);
   FEFaceValues<dim> fe_face_values_nbr (fe, quadrature,
                                         update_values);
   FESubfaceValues<dim> fe_subface_values (fe, quadrature,
                                           update_values | update_normal_vectors);
   FESubfaceValues<dim> fe_subface_values_nbr (fe, quadrature,
                                               update_values);
   
   QGauss<dim> q_cell (fe.degree + 1);
   FEValues<dim> fe_values (fe, q_cell, update_values);
   std::vector<double> sol_values(q_cell.size());
   
   std::vector<unsigned int> dof_indices_cell (fe_cell.dofs_per_cell);

   unsigned int n_q_points = quadrature.size();
   std::vector<double> face_values(n_q_points), face_values_nbr(n_q_points);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell0 = dof_handler_cell.begin_active();
   
   for(; cell != endc; ++cell, ++cell0)
   {
      cell0->get_dof_indices(dof_indices_cell);
      double& cell_shock_ind = shock_indicator (dof_indices_cell[0]);
      double& cell_jump_ind = jump_indicator (dof_indices_cell[0]);
      
      cell_shock_ind = 0;
      cell_jump_ind = 0;
      double inflow_measure = 0;
      
      // advection speed at cell center
      Point<dim> beta;
      advection_speed(cell->center(), beta);
      
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
         if (cell->at_boundary(f) == false)
         {
            if ((cell->neighbor(f)->level() == cell->level()) &&
                (cell->neighbor(f)->has_children() == false))
            {
               fe_face_values.reinit(cell, f);
               fe_face_values_nbr.reinit(cell->neighbor(f), cell->neighbor_of_neighbor(f));
               fe_face_values.get_function_values(solution, face_values);
               fe_face_values_nbr.get_function_values(solution, face_values_nbr);
               for(unsigned int q=0; q<n_q_points; ++q)
               {
                  int inflow_status = (beta * fe_face_values.normal_vector(q) < 0);
                  cell_shock_ind += inflow_status *
                                    (face_values[q] - face_values_nbr[q]) *
                                    fe_face_values.JxW(q);
                  cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                   fe_face_values.JxW(q);
                  inflow_measure += inflow_status * fe_face_values.JxW(q);
               }
               
            }
            else if ((cell->neighbor(f)->level() == cell->level()) &&
                     (cell->neighbor(f)->has_children() == true))
            {
               for (unsigned int subface=0; subface<cell->face(f)->n_children(); ++subface)
               {
                  fe_subface_values.reinit (cell, f, subface);
                  fe_face_values_nbr.reinit (cell->neighbor_child_on_subface (f, subface),
                                             cell->neighbor_of_neighbor(f));
                  fe_subface_values.get_function_values(solution, face_values);
                  fe_face_values_nbr.get_function_values(solution, face_values_nbr);
                  for(unsigned int q=0; q<n_q_points; ++q)
                  {
                     int inflow_status = (beta * fe_subface_values.normal_vector(q) < 0);
                     cell_shock_ind += inflow_status *
                                       (face_values[q] - face_values_nbr[q]) *
                                       fe_subface_values.JxW(q);
                     cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                      fe_face_values.JxW(q);
                     inflow_measure += inflow_status * fe_subface_values.JxW(q);
                  }
               }
            }
            else if (cell->neighbor_is_coarser(f))
            {
               fe_face_values.reinit(cell, f);
               fe_subface_values_nbr.reinit (cell->neighbor(f),
                                             cell->neighbor_of_coarser_neighbor(f).first,
                                             cell->neighbor_of_coarser_neighbor(f).second);
               fe_face_values.get_function_values(solution, face_values);
               fe_subface_values_nbr.get_function_values(solution, face_values_nbr);
               for(unsigned int q=0; q<n_q_points; ++q)
               {
                  int inflow_status = (beta * fe_face_values.normal_vector(q) < 0);
                  cell_shock_ind += inflow_status *
                                    (face_values[q] - face_values_nbr[q]) *
                                    fe_face_values.JxW(q);
                  cell_jump_ind += std::pow(face_values[q] - face_values_nbr[q], 2) *
                                   fe_face_values.JxW(q);
                  inflow_measure += inflow_status * fe_face_values.JxW(q);
               }
            }
         }
         else
         {
            // Boundary face
            // We dont do anything here since we assume solution is constant near
            // boundary.
         }
      
      // normalized shock indicator
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, sol_values);
      double cell_norm = 0;
      for(unsigned int q=0; q<q_cell.size(); ++q)
         cell_norm = std::max(cell_norm, std::fabs(sol_values[q]));
      
      double denominator = std::pow(cell->diameter(), 0.5*(fe.degree+1)) *
                           inflow_measure *
                           cell_norm;
      if(denominator > 1.0e-12)
      {
         cell_shock_ind = std::fabs(cell_shock_ind) / denominator;
         cell_shock_ind = (cell_shock_ind > 1.0) ? 1.0 : 0.0;
      }
      else
         cell_shock_ind = 0;
      
      dx = cell->diameter() / std::sqrt(2.0);
      cell_jump_ind = std::sqrt( cell_jump_ind / (4.0*dx) ) * cell->diameter();
   }
}
//------------------------------------------------------------------------------
// Slope limiter based on minmod
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::apply_limiter_TVD ()
{
   if(fe.degree == 0) return;

   QTrapez<1> q_trapez;
   QMidpoint<1> q_midpoint;
   QAnisotropic<dim> qrule_x (q_trapez, q_midpoint);
   FEValues<dim> fe_values_x (fe, qrule_x, update_values);
   QAnisotropic<dim> qrule_y (q_midpoint, q_trapez);
   FEValues<dim> fe_values_y (fe, qrule_y, update_values);

   std::vector<unsigned int> dof_indices     (fe.dofs_per_cell);
   std::vector<unsigned int> dof_indices_nbr (fe.dofs_per_cell);
   std::vector<double> face_values(2);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   double average_nbr, change, change_x, change_y;
   double db, df, DB, DF, DBx_new, DFx_new, DBy_new, DFy_new;

   for(unsigned int c=0; cell != endc; ++c, ++cell)
   {
      // mesh width: we assume cartesian mesh
      dx = cell->diameter() / std::sqrt(2.0);
      
      cell->get_dof_indices (dof_indices);

      // Limit x derivative
      fe_values_x.reinit(cell);
      fe_values_x.get_function_values(solution, face_values);

      df = db = 0;

      if(lcell[c] != endc)
      {
         //lcell[c]->get_dof_indices (dof_indices_nbr);
         //db = solution(dof_indices[0]) - solution(dof_indices_nbr[0]);
         average_nbr = compute_cell_average(lcell[c]);
         db = solution(dof_indices[0]) - average_nbr;
      }

      if(rcell[c] != endc)
      {
         //rcell[c]->get_dof_indices (dof_indices_nbr);
         //df = solution(dof_indices_nbr[0]) - solution(dof_indices[0]);
         average_nbr = compute_cell_average(rcell[c]);
         df = average_nbr - solution(dof_indices[0]);
      }

      DB = solution(dof_indices[0]) - face_values[0];
      DF = face_values[1] - solution(dof_indices[0]);

      DBx_new = minmod(DB, db, df);
      DFx_new = minmod(DF, db, df);
      change_x = std::fabs(DBx_new - DB) + std::fabs(DFx_new - DF);

      // Limit y derivative
      fe_values_y.reinit(cell);
      fe_values_y.get_function_values(solution, face_values);

      df = db = 0;

      if(bcell[c] != endc)
      {
         //bcell[c]->get_dof_indices (dof_indices_nbr);
         //db = solution(dof_indices[0]) - solution(dof_indices_nbr[0]);
         average_nbr = compute_cell_average(bcell[c]);
         db = solution(dof_indices[0]) - average_nbr;
      }

      if(tcell[c] != endc)
      {
         //tcell[c]->get_dof_indices (dof_indices_nbr);
         //df = solution(dof_indices_nbr[0]) - solution(dof_indices[0]);
         average_nbr = compute_cell_average(tcell[c]);
         df = average_nbr - solution(dof_indices[0]);
      }

      DB = solution(dof_indices[0]) - face_values[0];
      DF = face_values[1] - solution(dof_indices[0]);
      DBy_new = minmod(DB, db, df);
      DFy_new = minmod(DF, db, df);
      change_y = std::fabs(DBy_new - DB) + std::fabs(DFy_new - DF);

      // If limiter is active, then retain only linear part of solution
      change = 0.25*(change_x + change_y);

      if(change > 1.0e-10)
      {
         // Dont change cell average; zero other dofs
         for(unsigned int i=1; i<fe.dofs_per_cell; ++i)
            solution(dof_indices[i]) = 0;
         // NOTE: This part depends on the ordering of the basis functions.
         solution(dof_indices[1]) = 0.5*(DBx_new + DFx_new) / fe_values_x.shape_value(1,1);
         solution(dof_indices[fe.degree+1]) = 0.5*(DBy_new + DFy_new) 
                                              / fe_values_y.shape_value(fe.degree+1,1);
      }

   }
}

//------------------------------------------------------------------------------
// Limit the solution
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::apply_limiter ()
{
   if(limiter_type == tvd)
      apply_limiter_TVD ();
}

//------------------------------------------------------------------------------
// Compute min and max of cell average values
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::compute_min_max ()
{
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   sol_min =  1.0e20; 
   sol_max = -1.0e20;

   for(; cell != endc; ++cell)
   {
      cell->get_dof_indices (dof_indices);
      sol_min = std::min( sol_min, solution(dof_indices[0]) );
      sol_max = std::max( sol_max, solution(dof_indices[0]) );
   }
}

//------------------------------------------------------------------------------
// Solve the problem to convergence by RK time integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::solve ()
{
   RHSIntegrator<dim> rhs_integrator (dof_handler);
   setup_mesh_worker (rhs_integrator);
   compute_dt ();
   
   std::cout << "Solving by RK ...\n";

   unsigned int iter = 0;
   double time = 0;
   while (time < 2*M_PI && iter < 10000)
   {
      solution_old = solution;

      // 3-stage RK scheme
      for(unsigned int r=0; r<3; ++r)
      {
         assemble_rhs (rhs_integrator);
         
         for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
            solution(i) = a_rk[r] * solution_old(i) +
                          b_rk[r] * (solution(i) + dt * right_hand_side(i));

         apply_limiter ();
      }
      compute_min_max();
      
      ++iter; time += dt;

//      if(iter==1 || std::fmod(iter,10)==0)
//      {
//         compute_shock_indicator ();
//         refine_grid ();
//         compute_dt ();
//      }
      
      std::cout << "Iterations=" << iter
                << ", t = " << time 
                << ", min,max = " << sol_min << "  " << sol_max << endl;
      if(std::fmod(iter,50)==0) output_results(iter);
   }
   
}

//------------------------------------------------------------------------------
// Refine grid
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::refine_grid ()
{
   std::cout << "Refining grid ...\n";
   Vector<float> gradient_indicator (triangulation.n_active_cells());
   
   DerivativeApproximation::approximate_gradient (mapping,
                                                  dof_handler,
                                                  solution,
                                                  gradient_indicator);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      gradient_indicator(cell_no) *= std::pow(cell->diameter(), 1+1.0*dim/2);
   
   SolutionTransfer<dim, Vector<double> > soltrans(dof_handler);
   GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                      jump_indicator,
                                                      0.1, 0.5);
   
   unsigned int min_grid_level = 0;
   unsigned int max_grid_level = 5;
   if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
            cell = triangulation.begin_active(max_grid_level);
            cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
   
   for (typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active(min_grid_level);
        cell != triangulation.end_active(min_grid_level); ++cell)
      cell->clear_coarsen_flag ();
   
   // store solution on current mesh
   Vector<double> previous_solution;
   previous_solution = solution;
   
   triangulation.prepare_coarsening_and_refinement();
   soltrans.prepare_for_coarsening_and_refinement(previous_solution);
   triangulation.execute_coarsening_and_refinement ();
   
   setup_system ();
   
   // interpolate solution to new mesh
   soltrans.interpolate(previous_solution, solution);
   soltrans.clear();
}

//------------------------------------------------------------------------------
// Save results to file
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::output_results (const unsigned int cycle)
{
   // we copy average solution into "average"
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   std::vector<unsigned int> dof_indices_cell (fe_cell.dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      cell0= dof_handler_cell.begin_active(),
      endc = dof_handler.end();

   average = 0;
   
   for(; cell != endc; ++cell, ++cell0)
   {
      cell->get_dof_indices (dof_indices);
      cell0->get_dof_indices (dof_indices_cell);
      average(dof_indices_cell[0]) = solution(dof_indices[0]);
   }
   
   
   {
      // Output of the solution in
      std::string filename = "sol-" + Utilities::int_to_string(cycle,5) + ".vtk";
      std::cout << "Writing solution to <" << filename << ">" << std::endl;
      std::ofstream outfile (filename.c_str());
      
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "u");
      data_out.build_patches (fe.degree+1);
      data_out.write_vtk (outfile);
   }
   
   {
      // Output of the solution in
      std::string filename = "cell-" + Utilities::int_to_string(cycle,5) + ".vtk";
      std::cout << "Writing solution to <" << filename << ">" << std::endl;
      std::ofstream outfile (filename.c_str());
      
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler_cell);
      data_out.add_data_vector (average, "average");
      
      compute_shock_indicator ();
      data_out.add_data_vector (shock_indicator, "shock_indicator");
      data_out.add_data_vector (jump_indicator, "jump_indicator");
      
      data_out.build_patches (1);
      data_out.write_vtk (outfile);
   }
}

//------------------------------------------------------------------------------
// Actual computation starts from here
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::run ()
{
   GridGenerator::subdivided_hyper_cube (triangulation,40,-1.0,+1.0);
   setup_system ();
   set_initial_condition ();
   
   // Initial refinements
//   unsigned int n_refine_init = 1;
//   for(unsigned int i=0; i<n_refine_init; ++i)
//   {
//      compute_shock_indicator ();
//      refine_grid ();
//      set_initial_condition();
//   }
   output_results(0);
   
   solve ();
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main ()
{
   try
   {
      unsigned int degree = 1;
      LimiterType limiter_type = tvd;
      TestCase test_case = square;
      Mlim = 100.0;
      Step12<2> dgmethod(degree, limiter_type, test_case);
      dgmethod.run ();
   }
   catch (std::exception &exc)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
		          << std::endl;
      std::cerr << "Exception on processing: " << std::endl
		          << exc.what() << std::endl
		          << "Aborting!" << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   };
   
   return 0;
}


