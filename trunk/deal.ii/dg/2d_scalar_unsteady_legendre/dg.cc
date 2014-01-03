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
enum TestCase {expo, square};

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
            values[i] = std::pow(std::cos(2.0*M_PI*r), 2.0);
         else
            values[i] = 0.0;
      }
   else if(test_case == square)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         const Point<dim>& p = points[i];
         if(std::fabs(p(0)-0.5) < 0.25 && std::fabs(p(1)) < 0.25)
            values[i] = 1.0;
         else
            values[i] = 0.0;
      }
   else
      for (unsigned int i=0; i<values.size(); ++i)
         values[i] = 0.0;
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
      values[i]=0.0;
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
      void apply_limiter_old ();
      void apply_limiter ();
      void apply_limiter_TVD ();
      void refine_grid ();
      void compute_min_max ();
      void output_results (const unsigned int cycle);
      
      Triangulation<dim>   triangulation;
      const MappingQ1<dim> mapping;
   
      unsigned int         degree;
      FE_DGP<dim>          fe;
      DoFHandler<dim>      dof_handler;
      
      std::vector< Vector<double> > inv_mass_matrix;
      
      Vector<double>       solution;
      Vector<double>       solution_old;
      Vector<double>       average;
      Vector<double>       right_hand_side;
      double               dt;
      double               cfl;
      LimiterType          limiter_type;
      TestCase             test_case;
      double               sol_min, sol_max;
   
      std::vector<typename DoFHandler<dim>::active_cell_iterator>
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
      limiter_type (limiter_type),
      test_case (test_case)
{
   cfl = 0.4/(2.0*degree + 1.0);
}

//------------------------------------------------------------------------------
// Make dofs and allocate memory
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_system ()
{
   std::cout << "Allocating memory ...\n";

   dof_handler.distribute_dofs (fe);
   
   inv_mass_matrix.resize(triangulation.n_cells());
   for (unsigned int c=0; c<triangulation.n_cells(); ++c)
      inv_mass_matrix[c].reinit(fe.dofs_per_cell);

   solution.reinit (dof_handler.n_dofs());
   solution_old.reinit (dof_handler.n_dofs());
   average.reinit (dof_handler.n_dofs());
   right_hand_side.reinit (dof_handler.n_dofs());
   
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

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
         if (! cell->at_boundary(face_no))
         {
            const typename DoFHandler<dim>::cell_iterator
            neighbor = cell->neighbor(face_no);
            if (neighbor->active())
            {
               Point<dim> dr = neighbor->center() - cell->center();
               if(dr(0) < 0 && std::fabs(dr(1)) < EPS)
                  lcell[c] = neighbor;
               else if(dr(0) > 0 && std::fabs(dr(1)) < EPS)
                  rcell[c] = neighbor;
               else if(dr(1) < 0 && std::fabs(dr(0)) < EPS)
                  bcell[c] = neighbor;
               else if(dr(1) > 0 && std::fabs(dr(0)) < EPS)
                  tcell[c] = neighbor;
               else
               {
                  std::cout << "Did not find all neighbours\n";
                  exit(0);
               }
            }
            else
            {
               std::cout << "Not active cell\n";
                  exit(0);
            }
         }
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

   const unsigned int n_gauss_points = dof_handler.get_fe().degree+1;
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
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
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
      
      const double beta_n=beta * normals[point];
      if (beta_n>0)
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector(i) -= beta_n * 
                               sol[point] *
                               fe_v.shape_value(i,point) *
                               JxW[point];
      else
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector(i) -= beta_n *
                               g[point] *
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
      
      const double beta_n=beta * normals[point];
      if (beta_n>0)
      {
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= beta_n *
                                sol1[point] *
                                fe_v.shape_value(i,point) *
                                JxW[point];
         
         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += beta_n *
                                sol1[point] *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
      }
      else
      {
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= beta_n *
                                sol2[point] *
                                fe_v.shape_value(i,point) *
                                JxW[point];
         
         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += beta_n *
                                sol2[point] *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
      }
   }
}

//------------------------------------------------------------------------------
// Slope limiter based on minmod
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::apply_limiter_old ()
{
   if(fe.degree == 0) return;

   Assert (fe.degree == 1, ExcInternalError() );

   QMidpoint<dim> qrule;
   FEValues<dim> fe_values (fe, qrule, update_gradients);
   
   std::vector<unsigned int> dof_indices     (fe.dofs_per_cell);
   std::vector<unsigned int> dof_indices_nbr (fe.dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   // mesh width: we assume uniform cartesian mesh
   dx = cell->diameter() / std::sqrt(2.0);

   for(unsigned int c=0; cell != endc; ++c, ++cell)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices (dof_indices);
      const Tensor<1,dim>& grad_phi1 = fe_values.shape_grad (1,0);
      const Tensor<1,dim>& grad_phi2 = fe_values.shape_grad (2,0);
      Tensor<1,dim> grad_u = solution(dof_indices[1]) * grad_phi1
                           + solution(dof_indices[2]) * grad_phi2;

      double dal = 0, dar = 0;
      double dab = 0, dat = 0;

      if(lcell[c] != endc)
      {
         lcell[c]->get_dof_indices (dof_indices_nbr);
         dal = (solution(dof_indices[0]) - solution(dof_indices_nbr[0]))/dx;
      }

      if(rcell[c] != endc)
      {
         rcell[c]->get_dof_indices (dof_indices_nbr);
         dar = (solution(dof_indices_nbr[0]) - solution(dof_indices[0]))/dx;
      }

      if(bcell[c] != endc)
      {
         bcell[c]->get_dof_indices (dof_indices_nbr);
         dab = (solution(dof_indices[0]) - solution(dof_indices_nbr[0]))/dx;
      }

      if(tcell[c] != endc)
      {
         tcell[c]->get_dof_indices (dof_indices_nbr);
         dat = (solution(dof_indices_nbr[0]) - solution(dof_indices[0]))/dx;
      }

      double u_x = minmod(grad_u[0], dal, dar);
      double u_y = minmod(grad_u[1], dab, dat);
      solution(dof_indices[1]) = u_x / grad_phi1[0];
      solution(dof_indices[2]) = u_y / grad_phi2[1];
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

   // mesh width: we assume uniform cartesian mesh
   dx = cell->diameter() / std::sqrt(2.0);

   double change, change_x, change_y;
   double db, df, DB, DF, DBx_new, DFx_new, DBy_new, DFy_new;

   for(unsigned int c=0; cell != endc; ++c, ++cell)
   {
      cell->get_dof_indices (dof_indices);

      // Limit x derivative
      fe_values_x.reinit(cell);
      fe_values_x.get_function_values(solution, face_values);

      df = db = 0;

      if(lcell[c] != endc)
      {
         lcell[c]->get_dof_indices (dof_indices_nbr);
         db = solution(dof_indices[0]) - solution(dof_indices_nbr[0]);
      }

      if(rcell[c] != endc)
      {
         rcell[c]->get_dof_indices (dof_indices_nbr);
         df = solution(dof_indices_nbr[0]) - solution(dof_indices[0]);
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
         bcell[c]->get_dof_indices (dof_indices_nbr);
         db = solution(dof_indices[0]) - solution(dof_indices_nbr[0]);
      }

      if(tcell[c] != endc)
      {
         tcell[c]->get_dof_indices (dof_indices_nbr);
         df = solution(dof_indices_nbr[0]) - solution(dof_indices[0]);
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
      gradient_indicator(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);
   
   SolutionTransfer<dim, Vector<double> > soltrans(dof_handler);
   GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                      gradient_indicator,
                                                      0.3, 0.0);
   triangulation.prepare_coarsening_and_refinement();
   soltrans.prepare_for_coarsening_and_refinement(solution);
   triangulation.execute_coarsening_and_refinement ();
   dof_handler.distribute_dofs (fe);
   solution_old.reinit(dof_handler.n_dofs());
   soltrans.interpolate(solution, solution_old);
   soltrans.clear();
   solution.reinit(dof_handler.n_dofs());
   solution = solution_old;
}

//------------------------------------------------------------------------------
// Save results to file
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::output_results (const unsigned int cycle)
{
   // we copy average solution into "average"
   // This is inefficient but we dont care now.
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   average = 0;
   
   for(; cell != endc; ++cell)
   {
      cell->get_dof_indices (dof_indices);
      average(dof_indices[0]) = solution(dof_indices[0]);
   }
   
   // Output of the solution in
   std::string filename = "sol-" + Utilities::int_to_string(cycle,5) + ".vtk";
   std::cout << "Writing solution to <" << filename << ">" << std::endl;
   std::ofstream outfile (filename.c_str());
   
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "u");
   data_out.add_data_vector (average, "average");
   
   data_out.build_patches (fe.degree+1);
   
   data_out.write_vtk (outfile);
}

//------------------------------------------------------------------------------
// Actual computation starts from here
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::run ()
{
   GridGenerator::hyper_cube (triangulation,-1.0,+1.0);
   triangulation.refine_global (6);
   
   std::cout << "Number of active cells:       "
             << triangulation.n_active_cells()
             << std::endl;
   
   setup_system ();
   
   std::cout << "Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << std::endl;
   
   assemble_mass_matrix ();
   set_initial_condition ();
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
      Mlim = 0.0;
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


