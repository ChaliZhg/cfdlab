// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/vtk_io.h"

#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"

using namespace libMesh;

Real exact_solution(const Real x,
                    const Real y,
                    const Real t);

//------------------------------------------------------------------------------
Number exact_value(const Point& p,
                   const Parameters& parameters,
                   const std::string&,
                   const std::string&)
{
   return exact_solution(p(0), p(1), parameters.get<Real>("time"));
}

//------------------------------------------------------------------------------
void set_initial_condition(EquationSystems& es,
                           const std::string& system_name)
{
   libmesh_assert_equal_to (system_name, "Scalar-Convection");

   TransientExplicitSystem& system =
      es.get_system<TransientExplicitSystem>("Scalar-Convection");

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

   TransientExplicitSystem& system =
      es.get_system<TransientExplicitSystem> ("Scalar-Convection");

   FEType fe_type = system.variable_type (0);

   AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

   const DofMap& dof_map = system.get_dof_map ();

   DenseVector<Number> rhs;
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
         rhs.resize (n_dof);
      }

      system.rhs->add_vector (rhs, dof_indices);
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

   TransientExplicitSystem& system =
      equation_systems.add_system<TransientExplicitSystem> ("Scalar-Convection");
   system.add_variable ("u", FIRST, XYZ);
   system.attach_assemble_function (assemble);
   system.attach_init_function (set_initial_condition);

   equation_systems.init ();
   equation_systems.print_info ();

   VTKIO(mesh).write_equation_systems ("out_000.pvtu",
                                       equation_systems);

   return 0;
}
