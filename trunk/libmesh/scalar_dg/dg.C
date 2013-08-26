// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/vtk_io.h"

#include "libmesh/quadrature_gauss.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"

using namespace libMesh;

Real exact_solution(const Real x,
                    const Real y,
                    const Real t);

//------------------------------------------------------------------------------
Real exact_value(const Point& p,
                 const Parameters& parameters,
                 const std::string&,
                 const std::string&)
{
   return exact_solution(p(0), p(1), parameters.get<Real>("time"));
}

//------------------------------------------------------------------------------
Point advection_speed(const Point& p)
{
   Gradient beta;
   beta(0) = -p(1);
   beta(1) =  p(0);
   return beta;
}
//------------------------------------------------------------------------------
void set_initial_condition(EquationSystems& es,
                           const std::string& system_name)
{
   libmesh_assert_equal_to (system_name, "Scalar-Convection");

   TransientLinearImplicitSystem& system =
      es.get_system<TransientLinearImplicitSystem>("Scalar-Convection");

   es.parameters.set<Real> ("time") = system.time = 0;

   system.project_solution (exact_value, NULL, es.parameters);
}

//------------------------------------------------------------------------------
void assemble(EquationSystems& es,
              const std::string& system_name)
{
   libmesh_assert_equal_to (system_name, "Scalar-Convection");

   const MeshBase& mesh = es.get_mesh ();
   const unsigned int dim = mesh.mesh_dimension ();

   TransientLinearImplicitSystem& system =
      es.get_system<TransientLinearImplicitSystem> ("Scalar-Convection");

   const DofMap& dof_map = system.get_dof_map ();
   FEType fe_type = system.variable_type (0);
   AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

   QGauss qrule (dim, SECOND);
   fe->attach_quadrature_rule (&qrule);

   const std::vector<Real>& JxW = fe->get_JxW ();
   const std::vector<Point>& q_point = fe->get_xyz ();
   const std::vector<std::vector<Real> >& phi = fe->get_phi ();
   const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi ();

   // Reciprocal of time step
   const Real rdt = 1.0 / es.parameters.get<Real> ("dt");

   DenseMatrix<Real> Me;
   DenseVector<Real> Re;
   std::vector<dof_id_type> dof_indices;

   MeshBase::const_element_iterator el = mesh.active_local_elements_begin ();
   const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end ();

   for( ; el != end_el; ++el)
   {
      const Elem* elem = *el;
      dof_map.dof_indices (elem, dof_indices);
      fe->reinit (elem);
      {
         const unsigned int n_dof = dof_indices.size ();
         Re.resize (n_dof);
         Me.resize (n_dof, n_dof);
      }

      // Volume terms
      for(unsigned int q=0; q<qrule.n_points(); ++q)
      {
         for(unsigned int i=0; i<phi.size(); ++i)
            for(unsigned int j=0; j<phi.size(); ++j)
               Me(i,j) += rdt * phi[i][q] * phi[j][q] * JxW[q];

         Real sol = 0;
         for(unsigned int i=0; i<phi.size(); ++i)
            sol += system.current_solution(dof_indices[i]) * phi[i][q];

         Point beta = advection_speed (q_point[q]);
         for(unsigned int i=0; i<phi.size(); ++i)
            Re(i) += rdt * sol * phi[i][q] * JxW[q]
                      + sol * (beta * dphi[i][q]) * JxW[q];
      }

      // Face terms

      system.rhs->add_vector    (Re, dof_indices);
      system.matrix->add_matrix (Me, dof_indices);
   }

}

//------------------------------------------------------------------------------
int main (int argc, char** argv)
{
   LibMeshInit init(argc, argv);

   int dim = 2;

   Mesh mesh(init.comm());

   int ps = 10;
   Real halfwidth = 1.0;
   Real halfheight = 1.0;

   MeshTools::Generation::build_square(mesh,
                                       ps,
                                       ps,
                                      -halfwidth, halfwidth,
                                      -halfheight, halfheight,
                                       QUAD4);
   mesh.print_info();

   EquationSystems equation_systems (mesh);

   TransientLinearImplicitSystem& system =
      equation_systems.add_system<TransientLinearImplicitSystem> ("Scalar-Convection");
   system.add_variable ("u", FIRST, XYZ);
   system.attach_assemble_function (assemble);
   system.attach_init_function (set_initial_condition);

   equation_systems.init ();
   equation_systems.print_info ();

   VTKIO(mesh).write_equation_systems ("out_000.pvtu",
                                       equation_systems);

   const unsigned int n_time_steps = 100;
   const Real dt = 0.001;
   system.time = 0.0;

   for(unsigned int ts=0; ts<n_time_steps; ++ts)
   {
      system.time += dt;
      equation_systems.parameters.set<Real>("time") = system.time;
      equation_systems.parameters.set<Real>("dt") = dt;

      *system.old_local_solution = *system.current_local_solution;

      system.solve();

   }

   return 0;
}
