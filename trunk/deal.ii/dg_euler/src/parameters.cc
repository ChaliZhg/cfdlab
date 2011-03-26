#include "equation.h"
#include "parameters.h"

using namespace dealii;

   
namespace Parameters
{
   
   void Solver::declare_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("linear solver");
      {
         prm.declare_entry("output", "quiet",
                           Patterns::Selection("quiet|verbose"),
                           "State whether output from solver runs should be printed. "
                           "Choices are <quiet|verbose>.");
         prm.declare_entry("method", "gmres",
                           Patterns::Selection("gmres|direct"),
                           "The kind of solver for the linear system. "
                           "Choices are <gmres|direct>.");
         prm.declare_entry("residual", "1e-10",
                           Patterns::Double(),
                           "Linear solver residual");
         prm.declare_entry("max iters", "300",
                           Patterns::Integer(),
                           "Maximum solver iterations");
         prm.declare_entry("ilut fill", "2",
                           Patterns::Double(),
                           "Ilut preconditioner fill");
         prm.declare_entry("ilut absolute tolerance", "1e-9",
                           Patterns::Double(),
                           "Ilut preconditioner tolerance");
         prm.declare_entry("ilut relative tolerance", "1.1",
                           Patterns::Double(),
                           "Ilut relative tolerance");
         prm.declare_entry("ilut drop tolerance", "1e-10",
                           Patterns::Double(),
                           "Ilut drop tolerance");
      }
      prm.leave_subsection();
   }
   
   
   
   
   void Solver::parse_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("linear solver");
      {
         const std::string op = prm.get("output");
         if (op == "verbose")
            output = verbose;
         if (op == "quiet")
            output = quiet;
         
         const std::string sv = prm.get("method");
         if (sv == "direct")
            solver = direct;
         else if (sv == "gmres")
            solver = gmres;
         
         linear_residual = prm.get_double("residual");
         max_iterations  = prm.get_integer("max iters");
         ilut_fill       = prm.get_double("ilut fill");
         ilut_atol       = prm.get_double("ilut absolute tolerance");
         ilut_rtol       = prm.get_double("ilut relative tolerance");
         ilut_drop       = prm.get_double("ilut drop tolerance");
      }
      prm.leave_subsection();
   }
   
   
   
   void Refinement::declare_parameters (ParameterHandler &prm)
   {
      
      prm.enter_subsection("refinement");
      {
         prm.declare_entry("refinement", "true",
                           Patterns::Bool(),
                           "Whether to perform mesh refinement or not");
         prm.declare_entry("refinement fraction", "0.1",
                           Patterns::Double(),
                           "Fraction of high refinement");
         prm.declare_entry("unrefinement fraction", "0.1",
                           Patterns::Double(),
                           "Fraction of low unrefinement");
         prm.declare_entry("max elements", "1000000",
                           Patterns::Double(),
                           "maximum number of elements");
         prm.declare_entry("shock value", "4.0",
                           Patterns::Double(),
                           "value for shock indicator");
         prm.declare_entry("shock levels", "3.0",
                           Patterns::Double(),
                           "number of shock refinement levels");
      }
      prm.leave_subsection();
   }
   
   
   void Refinement::parse_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("refinement");
      {
         do_refine     = prm.get_bool ("refinement");
         shock_val     = prm.get_double("shock value");
         shock_levels  = prm.get_double("shock levels");
      }
      prm.leave_subsection();
   }
   
   
   
   void Flux::declare_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("flux");
      {
         prm.declare_entry("flux", "lxf",
                           Patterns::Selection("lxf|sw|kfvs"),
                           "Numerical flux: lxf | sw | kfvs");
         prm.declare_entry("stab", "mesh",
                           Patterns::Selection("constant|mesh"),
                           "Whether to use a constant stabilization parameter or "
                           "a mesh-dependent one");
         prm.declare_entry("stab value", "1",
                           Patterns::Double(),
                           "alpha stabilization");
      }
      prm.leave_subsection();
   }
   
   
   void Flux::parse_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("flux");
      {
         const std::string flux = prm.get("flux");
         if(flux == "lxf")
            flux_type = lxf;
         else if(flux == "sw")
            flux_type = sw;
         else if(flux == "kfvs")
            flux_type = kfvs;
         else
            AssertThrow (false, ExcNotImplemented());

         const std::string stab = prm.get("stab");
         if (stab == "constant")
            stabilization_kind = constant;
         else if (stab == "mesh")
            stabilization_kind = mesh_dependent;
         else
            AssertThrow (false, ExcNotImplemented());
         
         stabilization_value = prm.get_double("stab value");
      }
      prm.leave_subsection();
   }
   
   
   
   void Output::declare_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("output");
      {
         prm.declare_entry("schlieren plot", "false",
                           Patterns::Bool (),
                           "Whether or not to produce schlieren plots");
         prm.declare_entry("step", "-1",
                           Patterns::Double(),
                           "Output once per this period");
         prm.declare_entry("format", "vtk",
                           Patterns::Selection("vtk|tecplot"),
                           "Output format for visualization: vtk, tecplot");
      }
      prm.leave_subsection();
   }
   
   
   
   void Output::parse_parameters (ParameterHandler &prm)
   {
      prm.enter_subsection("output");
      {
         schlieren_plot = prm.get_bool("schlieren plot");
         output_step = prm.get_double("step");
         output_format = prm.get("format");
      }
      prm.leave_subsection();
   }
   
   
   
   template <int dim>
   AllParameters<dim>::BoundaryConditions::BoundaryConditions ()
   :
   values (EulerEquations<dim>::n_components)
   {}
   
   
   template <int dim>
   AllParameters<dim>::AllParameters ()
   :
   initial_conditions (EulerEquations<dim>::n_components)
   {}
   
   
   template <int dim>
   void
   AllParameters<dim>::declare_parameters (ParameterHandler &prm)
   {
      prm.declare_entry("mesh type", "gmsh",
                        Patterns::Selection("ucd|gmsh"),
                        "Mesh file type");
      
      prm.declare_entry("mesh file", "grid.msh",
                        Patterns::Anything(),
                        "Mesh file name");
      
      prm.declare_entry("diffusion power", "2.0",
                        Patterns::Double(),
                        "power of mesh size for diffusion");
      
      prm.enter_subsection("time stepping");
      {
         prm.declare_entry("time step", "0.1",
                           Patterns::Double(0),
                           "simulation time step");
         prm.declare_entry("final time", "10.0",
                           Patterns::Double(0),
                           "simulation end time");
         prm.declare_entry("theta scheme value", "0.5",
                           Patterns::Double(0,1),
                           "value for theta that interpolated between explicit "
                           "Euler (theta=0), Crank-Nicolson (theta=0.5), and "
                           "implicit Euler (theta=1).");
         prm.declare_entry("nonlinear iterations", "10",
                           Patterns::Integer(),
                           "maximum non-linear iterations");
      }
      prm.leave_subsection();
      
      
      for (unsigned int b=0; b<max_n_boundaries; ++b)
      {
         prm.enter_subsection("boundary_" +
                              Utilities::int_to_string(b));
         {
            prm.declare_entry("type", "outflow",
                              Patterns::Selection("slip|inflow|outflow|pressure|farfield"),
                              "<slip|inflow|outflow|pressure|farfield>");
            
            for (unsigned int di=0; di<EulerEquations<dim>::n_components; ++di)
            {
               prm.declare_entry("w_" + Utilities::int_to_string(di) +
                                 " value", "0.0",
                                 Patterns::Anything(),
                                 "expression in x,y,z");
            }
         }
         prm.leave_subsection();
      }
      
      prm.enter_subsection("initial condition");
      {
         for (unsigned int di=0; di<EulerEquations<dim>::n_components; ++di)
            prm.declare_entry("w_" + Utilities::int_to_string(di) + " value",
                              "0.0",
                              Patterns::Anything(),
                              "expression in x,y,z");
      }
      prm.leave_subsection();
      
      Parameters::Solver::declare_parameters (prm);
      Parameters::Refinement::declare_parameters (prm);
      Parameters::Flux::declare_parameters (prm);
      Parameters::Output::declare_parameters (prm);
   }
   
   
   template <int dim>
   void
   AllParameters<dim>::parse_parameters (ParameterHandler &prm)
   {
      mesh_type = prm.get("mesh type");
      mesh_filename = prm.get("mesh file");
      diffusion_power = prm.get_double("diffusion power");
      
      prm.enter_subsection("time stepping");
      {
         time_step = prm.get_double("time step");
         if (time_step == 0)
         {
            is_stationary = true;
            time_step = 1.0;
            final_time = 1.0e20;
         }
         else
            is_stationary = false;
         
         final_time = prm.get_double("final time");
         theta = prm.get_double("theta scheme value");
         max_nonlin_iter = prm.get_integer("nonlinear iterations");
      }
      prm.leave_subsection();
      
      for (unsigned int boundary_id=0; boundary_id<max_n_boundaries;
           ++boundary_id)
      {
         prm.enter_subsection("boundary_" +
                              Utilities::int_to_string(boundary_id));
         {
            std::vector<std::string> expressions(EulerEquations<dim>::n_components, "0.0");
            
            std::string boundary_type = prm.get("type");

            if (boundary_type == "slip")
               boundary_conditions[boundary_id].kind
               = EulerEquations<dim>::no_penetration_boundary;
            else if (boundary_type == "inflow")
               boundary_conditions[boundary_id].kind
               = EulerEquations<dim>::inflow_boundary;
            else if (boundary_type == "pressure")
               boundary_conditions[boundary_id].kind
               = EulerEquations<dim>::pressure_boundary;
            else if (boundary_type == "outflow")
               boundary_conditions[boundary_id].kind
               = EulerEquations<dim>::outflow_boundary;
            else if (boundary_type == "farfield")
               boundary_conditions[boundary_id].kind
               = EulerEquations<dim>::farfield_boundary;
            else
               AssertThrow (false, ExcNotImplemented());
            
            for (unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
            {               
               expressions[c] = prm.get("w_" + Utilities::int_to_string(c) +
                                         " value");
            }
            
            boundary_conditions[boundary_id].values.initialize 
               (FunctionParser<dim>::default_variable_names(),
                expressions,
                std::map<std::string, double>());
         }
         prm.leave_subsection();
      }
      
      prm.enter_subsection("initial condition");
      {
         std::vector<std::string> expressions (EulerEquations<dim>::n_components,
                                               "0.0");
         for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++)
            expressions[di] = prm.get("w_" + Utilities::int_to_string(di) +
                                      " value");
         initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                        expressions,
                                        std::map<std::string, double>());
      }
      prm.leave_subsection();
      
      Parameters::Solver::parse_parameters (prm);
      Parameters::Refinement::parse_parameters (prm);
      Parameters::Flux::parse_parameters (prm);
      Parameters::Output::parse_parameters (prm);
   }
}

// To avoid linking errors
template struct Parameters::AllParameters<2>;
