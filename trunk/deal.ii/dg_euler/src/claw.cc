#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>

#include <lac/vector.h>
#include <lac/compressed_sparsity_pattern.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_in.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>

#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>

#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <numerics/solution_transfer.h>

#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_precondition.h>
#include <lac/trilinos_solver.h>

#include <lac/solver_gmres.h>
#include <lac/precondition_block.h>

#include <Sacado.hpp>


#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "claw.h"

using namespace dealii;

// @sect4{ConservationLaw::ConservationLaw}
//
// There is nothing much to say about
// the constructor. Essentially, it
// reads the input file and fills the
// parameter object with the parsed
// values:
template <int dim>
ConservationLaw<dim>::ConservationLaw (const char *input_filename)
   :
   mapping (),
   fe (FE_DGQ<dim>(1), EulerEquations<dim>::n_components),
   dof_handler (triangulation),
   verbose_cout (std::cout, false)
{
   ParameterHandler prm;
   Parameters::AllParameters<dim>::declare_parameters (prm);
   
   prm.read_input (input_filename);
   parameters.parse_parameters (prm);
   
   verbose_cout.set_condition (parameters.output == Parameters::Solver::verbose);
}



// @sect4{ConservationLaw::setup_system}
//
// The following (easy) function is called
// each time the mesh is changed. All it
// does is to resize the Trilinos matrix
// according to a sparsity pattern that we
// generate as in all the previous tutorial
// programs.
template <int dim>
void ConservationLaw<dim>::setup_system ()
{
   DoFRenumbering::Cuthill_McKee (dof_handler);

   CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
   DoFTools::make_flux_sparsity_pattern (dof_handler, c_sparsity);

   SparsityPattern      sparsity_pattern;
   sparsity_pattern.copy_from(c_sparsity);

   // Visualize sparsity pattern
   //std::ofstream out ("sparsity_pattern.1");
   //sparsity_pattern.print_gnuplot (out);
   //abort ();
   
   system_matrix.reinit (sparsity_pattern);
}

//------------------------------------------------------------------------------
// Create mesh worker for integration
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::setup_mesh_worker (Integrator<dim>& integrator)
{   
   const unsigned int n_gauss_points = fe.degree + 1;
   integrator.info_box.initialize_gauss_quadrature(n_gauss_points,
                                                   n_gauss_points,
                                                   n_gauss_points);
   
   integrator.info_box.initialize_update_flags ();
   integrator.info_box.add_update_flags_all (update_values | 
                                             update_quadrature_points |
                                             update_JxW_values);
   integrator.info_box.add_update_flags_cell     (update_gradients);
   integrator.info_box.add_update_flags_boundary (update_normal_vectors);
   integrator.info_box.add_update_flags_face     (update_normal_vectors);
   
   integrator.info_box.initialize (fe, mapping);
   
   integrator.assembler.initialize (system_matrix, right_hand_side);
}

//------------------------------------------------------------------------------
// Contribution of volume integral terms
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_cell_term (DoFInfo& dinfo, 
                                           CellInfo& info)
{
   FullMatrix<double>& local_matrix = dinfo.matrix(0).matrix;
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;

   const FEValuesBase<dim>& fe_v    = info.fe_values();
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   const unsigned int n_q_points    = fe_v.n_quadrature_points;
   
   Table<2,Sacado::Fad::DFad<double> >
   W (n_q_points, EulerEquations<dim>::n_components);
   
   Table<2,double>
   W_old (n_q_points, EulerEquations<dim>::n_components);
   
   Table<2,Sacado::Fad::DFad<double> >
   W_theta (n_q_points, EulerEquations<dim>::n_components);
   
   Table<3,Sacado::Fad::DFad<double> >
   grad_W (n_q_points, EulerEquations<dim>::n_components, dim);
      
   // Next, we have to define the independent
   // variables that we will try to determine
   // by solving a Newton step. These
   // independent variables are the values of
   // the local degrees of freedom which we
   // extract here:
   std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
   
   // The next step incorporates all the
   // magic: we declare a subset of the
   // autodifferentiation variables as
   // independent degrees of freedom, whereas
   // all the other ones remain dependent
   // functions. These are precisely the local
   // degrees of freedom just extracted. All
   // calculations that reference them (either
   // directly or indirectly) will accumulate
   // sensitivies with respect to these
   // variables.
   //
   // In order to mark the variables as
   // independent, the following does the
   // trick, marking
   // <code>independent_local_dof_values[i]</code>
   // as the $i$th independent variable out of
   // a total of <code>dofs_per_cell</code>:
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i].diff (i, dofs_per_cell);
   
   // After all these declarations, let us
   // actually compute something. First, the
   // values of <code>W</code>,
   // <code>W_old</code>,
   // <code>W_theta</code>, and
   // <code>grad_W</code>, which we can
   // compute from the local DoF values by
   // using the formula $W(x_q)=\sum_i \mathbf
   // W_i \Phi_i(x_q)$, where $\mathbf W_i$ is
   // the $i$th entry of the (local part of
   // the) solution vector, and $\Phi_i(x_q)$
   // the value of the $i$th vector-valued
   // shape function evaluated at quadrature
   // point $x_q$. The gradient can be
   // computed in a similar way.
   //
   // Ideally, we could compute this
   // information using a call into something
   // like FEValues::get_function_values and
   // FEValues::get_function_grads, but since
   // (i) we would have to extend the FEValues
   // class for this, and (ii) we don't want
   // to make the entire
   // <code>old_solution</code> vector fad
   // types, only the local cell variables, we
   // explicitly code the loop above. Before
   // this, we add another loop that
   // initializes all the fad variables to
   // zero:
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
      {
         W[q][c]       = 0;
         W_old[q][c]   = 0;
         W_theta[q][c] = 0;
         for (unsigned int d=0; d<dim; ++d)
            grad_W[q][c][d] = 0;
      }
   
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;
         
         W[q][c] += independent_local_dof_values[i] *
                    fe_v.shape_value_component(i, q, c);
         W_old[q][c] += old_solution(dof_indices[i]) *
                        fe_v.shape_value_component(i, q, c);
         W_theta[q][c] += (parameters.theta *
                           independent_local_dof_values[i]
                           +
                           (1-parameters.theta) *
                           old_solution(dof_indices[i])) *
                           fe_v.shape_value_component(i, q, c);
         
         for (unsigned int d = 0; d < dim; d++)
            grad_W[q][c][d] += independent_local_dof_values[i] *
                               fe_v.shape_grad_component(i, q, c)[d];
      }
   
   
   // Next, in order to compute the cell
   // contributions, we need to evaluate
   // $F(\tilde{\mathbf w})$ and
   // $G(\tilde{\mathbf w})$ at all quadrature
   // points. To store these, we also need to
   // allocate a bit of memory. Note that we
   // compute the flux matrices and right hand
   // sides in terms of autodifferentiation
   // variables, so that the Jacobian
   // contributions can later easily be
   // computed from it:
   typedef Sacado::Fad::DFad<double> FluxMatrix[EulerEquations<dim>::n_components][dim];
   FluxMatrix *flux = new FluxMatrix[n_q_points];
   
   typedef Sacado::Fad::DFad<double> ForcingVector[EulerEquations<dim>::n_components];
   ForcingVector *forcing = new ForcingVector[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      EulerEquations<dim>::compute_flux_matrix (W_theta[q], flux[q]);
      EulerEquations<dim>::compute_forcing_vector (W_theta[q], forcing[q]);
   }
   
   
   // We now have all of the pieces in place,
   // so perform the assembly.  We have an
   // outer loop through the components of the
   // system, and an inner loop over the
   // quadrature points, where we accumulate
   // contributions to the $i$th residual
   // $F_i$. The general formula for this
   // residual is given in the introduction
   // and at the top of this function. We can,
   // however, simplify it a bit taking into
   // account that the $i$th (vector-valued)
   // test function $\mathbf{z}_i$ has in
   // reality only a single nonzero component
   // (more on this topic can be found in the
   // @ref vector_valued module). It will be
   // represented by the variable
   // <code>component_i</code> below. With
   // this, the residual term can be
   // re-written as $F_i =
   // \left(\frac{(\mathbf{w}_{n+1} -
   // \mathbf{w}_n)_{\text{component\_i}}}{\delta
   // t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K$
   // $- \sum_{d=1}^{\text{dim}}
   // \left(\mathbf{F}
   // (\tilde{\mathbf{w}})_{\text{component\_i},d},
   // \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}}
   // {\partial x_d}\right)_K$ $+
   // \sum_{d=1}^{\text{dim}} h^{\eta}
   // \left(\frac{\partial
   // \mathbf{w}_{\text{component\_i}}}{\partial
   // x_d} , \frac{\partial
   // (\mathbf{z}_i)_{\text{component\_i}}}{\partial
   // x_d} \right)_K$
   // $-(\mathbf{G}(\tilde{\mathbf{w}}
   // )_{\text{component\_i}},
   // (\mathbf{z}_i)_{\text{component\_i}})_K$,
   // where integrals are understood to be
   // evaluated through summation over
   // quadrature points.
   //
   // We initialy sum all contributions of the
   // residual in the positive sense, so that
   // we don't need to negative the Jacobian
   // entries.  Then, when we sum into the
   // <code>right_hand_side</code> vector,
   // we negate this residual.
   for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
   {
      Sacado::Fad::DFad<double> F_i = 0;
      
      const unsigned int
      component_i = fe_v.get_fe().system_to_component_index(i).first;
      
      // The residual for each row (i) will be accumulating
      // into this fad variable.  At the end of the assembly
      // for this row, we will query for the sensitivities
      // to this variable and add them into the Jacobian.
      
      for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
      {
         if (parameters.is_stationary == false)
            F_i += 1.0 / parameters.time_step *
                   (W[point][component_i] - W_old[point][component_i]) *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
         
         for (unsigned int d=0; d<dim; d++)
            F_i -= flux[point][component_i][d] *
                   fe_v.shape_grad_component(i, point, component_i)[d] *
                   fe_v.JxW(point);
         
         // Diffusion term for shocks
         if(parameters.diffusion_power > 0.0)
            for (unsigned int d=0; d<dim; d++)
               F_i += 1.0*std::pow(fe_v.get_cell()->diameter(),
                                 parameters.diffusion_power) *
                     grad_W[point][component_i][d] *
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.JxW(point);
         
         F_i -= forcing[point][component_i] *
                fe_v.shape_value_component(i, point, component_i) *
                fe_v.JxW(point);
      }
      
      // At the end of the loop, we have to
      // add the sensitivities to the
      // matrix and subtract the residual
      // from the right hand side. Trilinos
      // FAD data type gives us access to
      // the derivatives using
      // <code>F_i.fastAccessDx(k)</code>,
      // so we store the data in a
      // temporary array. This information
      // about the whole row of local dofs
      // is then added to the Trilinos
      // matrix at once (which supports the
      // data types we have chosen).
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
         local_matrix (i, k) += F_i.fastAccessDx(k);
      }
      local_vector (i) -= F_i.val();
   }
   
   delete[] forcing;
   delete[] flux;
   
}


//------------------------------------------------------------------------------
// Contribution from boundary faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_boundary_term (DoFInfo& dinfo, 
                                                    CellInfo& info)
{
   FullMatrix<double>& local_matrix = dinfo.matrix(0).matrix;
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo.indices;
   const unsigned int& face_no = dinfo.face_number;
   const double& face_diameter = dinfo.face->diameter();
   const unsigned int& boundary_id = dinfo.face->boundary_indicator();
   
   const FEValuesBase<dim>& fe_v = info.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   std::vector<Sacado::Fad::DFad<double> >
   independent_local_dof_values (dofs_per_cell);
   
   const unsigned int n_independent_variables = dofs_per_cell;
   
   for (unsigned int i = 0; i < dofs_per_cell; i++)
   {
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
      independent_local_dof_values[i].diff(i, n_independent_variables);
   }
   
   // Next, we need to define the values of
   // the conservative variables $\tilde
   // {\mathbf W}$ on this side of the face
   // ($\tilde {\mathbf W}^+$) and on the
   // opposite side ($\tilde {\mathbf
   // W}^-$). The former can be computed in
   // exactly the same way as in the previous
   // function, but note that the
   // <code>fe_v</code> variable now is of
   // type FEFaceValues or FESubfaceValues:
   Table<2,Sacado::Fad::DFad<double> >
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);
   
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][component_i] += (parameters.theta *
                                   independent_local_dof_values[i]
                                   +
                                   (1.0-parameters.theta) *
                                   old_solution(dof_indices[i])) *
                                   fe_v.shape_value_component(i, q, component_i);
      }
   

   // On the other hand, if this is an
   // external boundary face, then the values
   // of $W^-$ will be either functions of
   // $W^+$, or they will be prescribed,
   // depending on the kind of boundary
   // condition imposed here.
   //
   // To start the evaluation, let us ensure
   // that the boundary id specified for this
   // boundary is one for which we actually
   // have data in the parameters
   // object. Next, we evaluate the function
   // object for the inhomogeneity.  This is a
   // bit tricky: a given boundary might have
   // both prescribed and implicit values.  If
   // a particular component is not
   // prescribed, the values evaluate to zero
   // and are ignored below.
   //
   // The rest is done by a function that
   // actually knows the specifics of Euler
   // equation boundary conditions. Note that
   // since we are using fad variables here,
   // sensitivities will be updated
   // appropriately, a process that would
   // otherwise be tremendously complicated.

   Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
           ExcIndexRange (boundary_id, 0,
                          Parameters::AllParameters<dim>::max_n_boundaries));
   
   std::vector<Vector<double> >
   boundary_values(n_q_points, Vector<double>(EulerEquations<dim>::n_components));
   parameters.boundary_conditions[boundary_id]
   .values.vector_value_list(fe_v.get_quadrature_points(),
                             boundary_values);
   
   
   typename EulerEquations<dim>::BoundaryKind boundary_kind = 
      parameters.boundary_conditions[boundary_id].kind;

   for (unsigned int q = 0; q < n_q_points; q++)
      EulerEquations<dim>::compute_Wminus (boundary_kind,
                                           fe_v.normal_vector(q),
                                           Wplus[q],
                                           boundary_values[q],
                                           Wminus[q]);
   
   
   // Now that we have $\mathbf w^+$ and
   // $\mathbf w^-$, we can go about computing
   // the numerical flux function $\mathbf
   // H(\mathbf w^+,\mathbf w^-, \mathbf n)$
   // for each quadrature point. Before
   // calling the function that does so, we
   // also need to determine the
   // Lax-Friedrich's stability parameter:
   typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   if(boundary_kind == EulerEquations<dim>::farfield_boundary)
   {
      for (unsigned int q=0; q<n_q_points; ++q)
         EulerEquations<dim>::steger_warming_flux(fe_v.normal_vector(q),
                                                  Wplus[q], Wminus[q],
                                                  normal_fluxes[q]);
   }
   /*else if(boundary_kind == EulerEquations<dim>::no_penetration_boundary)
   {
      for (unsigned int q=0; q<n_q_points; ++q)
         EulerEquations<dim>::no_penetration_flux(fe_v.normal_vector(q),
                                                  Wminus[q],
                                                  normal_fluxes[q]);
   }*/
   else
   {
      for (unsigned int q=0; q<n_q_points; ++q)
         numerical_normal_flux(fe_v.normal_vector(q),
                               Wplus[q], 
                               Wminus[q],
                               normal_fluxes[q]);
   }
   
   // Now assemble the face term in exactly
   // the same way as for the cell
   // contributions in the previous
   // function. The only difference is that if
   // this is an internal face, we also have
   // to take into account the sensitivies of
   // the residual contributions to the
   // degrees of freedom on the neighboring
   // cell:
   for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
         }
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix (i,k) += F_i.fastAccessDx(k);
         
         local_vector (i) -= F_i.val();
      }
   
   delete[] normal_fluxes;   
}



//------------------------------------------------------------------------------
// Contribution from interior faces
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                                CellInfo& info1, CellInfo& info2)
{
   FullMatrix<double>& local_matrix11 = dinfo1.matrix(0,false).matrix;
   FullMatrix<double>& local_matrix12 = dinfo1.matrix(0,true).matrix;
   Vector<double>& local_vector   = dinfo1.vector(0).block(0);
   std::vector<unsigned int>& dof_indices = dinfo1.indices;
   const unsigned int& face_no = dinfo1.face_number;
   const double& face_diameter = dinfo1.face->diameter();
   
   FullMatrix<double>& local_matrix_neighbor22 = dinfo2.matrix(0,false).matrix;
   FullMatrix<double>& local_matrix_neighbor21 = dinfo2.matrix(0,true).matrix;
   Vector<double>& local_vector_neighbor   = dinfo2.vector(0).block(0);
   std::vector<unsigned int>& dof_indices_neighbor = dinfo2.indices;
   const unsigned int& face_no_neighbor = dinfo2.face_number;

   const FEValuesBase<dim>& fe_v = info1.fe_values();
   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
   
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();
   const unsigned int dofs_per_cell_neighbor = fe_v_neighbor.dofs_per_cell;

   std::vector<Sacado::Fad::DFad<double> >
   independent_local_dof_values (dofs_per_cell),
   independent_neighbor_dof_values (dofs_per_cell_neighbor);
   
   const unsigned int n_independent_variables = dofs_per_cell + 
                                                dofs_per_cell_neighbor;
   
   for (unsigned int i = 0; i < dofs_per_cell; i++)
   {
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
      independent_local_dof_values[i].diff(i, n_independent_variables);
   }
   
   for (unsigned int i = 0; i < dofs_per_cell_neighbor; i++)
   {
      independent_neighbor_dof_values[i] = current_solution(dof_indices_neighbor[i]);
      independent_neighbor_dof_values[i].diff(i+dofs_per_cell, n_independent_variables);
   }
   
   
   // Next, we need to define the values of
   // the conservative variables $\tilde
   // {\mathbf W}$ on this side of the face
   // ($\tilde {\mathbf W}^+$) and on the
   // opposite side ($\tilde {\mathbf
   // W}^-$). The former can be computed in
   // exactly the same way as in the previous
   // function, but note that the
   // <code>fe_v</code> variable now is of
   // type FEFaceValues or FESubfaceValues:
   Table<2,Sacado::Fad::DFad<double> >
   Wplus (n_q_points, EulerEquations<dim>::n_components),
   Wminus (n_q_points, EulerEquations<dim>::n_components);
   
   for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
         Wplus[q][component_i] += (parameters.theta *
                                   independent_local_dof_values[i]
                                   +
                                   (1.0-parameters.theta) *
                                   old_solution(dof_indices[i])) *
                                   fe_v.shape_value_component(i, q, component_i);
      }
   
   // Computing $\tilde {\mathbf W}^-$ is a
   // bit more complicated. If this is an
   // internal face, we can compute it as
   // above by simply using the independent
   // variables from the neighbor:
      for (unsigned int q=0; q<n_q_points; ++q)
         for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
         {
            const unsigned int component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            Wminus[q][component_i] += (parameters.theta *
                                       independent_neighbor_dof_values[i]
                                       +
                                       (1.0-parameters.theta) *
                                       old_solution(dof_indices_neighbor[i]))*
                                       fe_v_neighbor.shape_value_component(i, q, component_i);
         }
   
   
   // Now that we have $\mathbf w^+$ and
   // $\mathbf w^-$, we can go about computing
   // the numerical flux function $\mathbf
   // H(\mathbf w^+,\mathbf w^-, \mathbf n)$
   // for each quadrature point. Before
   // calling the function that does so, we
   // also need to determine the
   // Lax-Friedrich's stability parameter:
   typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
   NormalFlux *normal_fluxes = new NormalFlux[n_q_points];
   
   for (unsigned int q=0; q<n_q_points; ++q)
      numerical_normal_flux(fe_v.normal_vector(q),
                            Wplus[q], 
                            Wminus[q],
                            normal_fluxes[q]);
   
   // Now assemble the face term in exactly
   // the same way as for the cell
   // contributions in the previous
   // function. The only difference is that if
   // this is an internal face, we also have
   // to take into account the sensitivies of
   // the residual contributions to the
   // degrees of freedom on the neighboring
   // cell:
   for (unsigned int i=0; i<dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v.get_fe().system_to_component_index(i).first;
            
            F_i += normal_fluxes[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
         }
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix11 (i,k) += F_i.fastAccessDx(k);

            for (unsigned int k=0; k<dofs_per_cell_neighbor; ++k)
               local_matrix12(i,k) += F_i.fastAccessDx(dofs_per_cell+k);
         
         local_vector(i) -= F_i.val();
      }
   
   // Contributions to neighbouring cell
   for (unsigned int i=0; i<dofs_per_cell_neighbor; ++i)
      if (fe_v_neighbor.get_fe().has_support_on_face(i, face_no_neighbor) == true)
      {
         Sacado::Fad::DFad<double> F_i = 0;
         
         for (unsigned int point=0; point<n_q_points; ++point)
         {
            const unsigned int
            component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            
            F_i -= normal_fluxes[point][component_i] *
                   fe_v_neighbor.shape_value_component(i, point, component_i) *
                   fe_v_neighbor.JxW(point);
         }
         
         for (unsigned int k=0; k<dofs_per_cell_neighbor; ++k)
            local_matrix_neighbor22 (i,k) += F_i.fastAccessDx(dofs_per_cell+k);         
         
         for (unsigned int k=0; k<dofs_per_cell; ++k)
            local_matrix_neighbor21 (i,k) += F_i.fastAccessDx(k);
         
         local_vector_neighbor (i) -= F_i.val();
      }
   
   delete[] normal_fluxes;
}

template <int dim>
void ConservationLaw<dim>::assemble_system (Integrator<dim>& integrator)
{
   MeshWorker::integration_loop<dim, dim>
   (dof_handler.begin_active(),
    dof_handler.end(),
    integrator.dof_info, 
    integrator.info_box,
    boost::bind(&ConservationLaw<dim>::integrate_cell_term, 
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_boundary_term,
                this, _1, _2),
    boost::bind(&ConservationLaw<dim>::integrate_face_term,
                this, _1, _2, _3, _4),
    integrator.assembler, true);
   
   system_matrix.compress ();
}


// @sect4{ConservationLaw::solve}
//
// Here, we actually solve the linear system,
// using either of Trilinos' Aztec or Amesos
// linear solvers. The result of the
// computation will be written into the
// argument vector passed to this
// function. The result is a pair of number
// of iterations and the final linear
// residual.

template <int dim>
std::pair<unsigned int, double>
ConservationLaw<dim>::solve (Vector<double> &newton_update)
{
   switch (parameters.solver)
   {
         // If the parameter file specified
         // that a direct solver shall be
         // used, then we'll get here. The
         // process is straightforward, since
         // deal.II provides a wrapper class
         // to the Amesos direct solver within
         // Trilinos. All we have to do is to
         // create a solver control object
         // (which is just a dummy object
         // here, since we won't perform any
         // iterations), and then create the
         // direct solver object. When
         // actually doing the solve, note
         // that we don't pass a
         // preconditioner. That wouldn't make
         // much sense for a direct solver
         // anyway.  At the end we return the
         // solver control statistics &mdash;
         // which will tell that no iterations
         // have been performed and that the
         // final linear residual is zero,
         // absent any better information that
         // may be provided here:
      case Parameters::Solver::direct:
      {
         SolverControl solver_control (1,0);
         TrilinosWrappers::SolverDirect direct (solver_control,
                                                parameters.output ==
                                                Parameters::Solver::verbose);
         
         direct.solve (system_matrix, newton_update, right_hand_side);
         
         return std::pair<unsigned int, double> (solver_control.last_step(),
                                                 solver_control.last_value());
      }
         
         // Likewise, if we are to use an
         // iterative solver, we use Aztec's
         // GMRES solver. We could use the
         // Trilinos wrapper classes for
         // iterative solvers and
         // preconditioners here as well, but
         // we choose to use an Aztec solver
         // directly. For the given problem,
         // Aztec's internal preconditioner
         // implementations are superior over
         // the ones deal.II has wrapper
         // classes to, so we use ILU-T
         // preconditioning within the AztecOO
         // solver and set a bunch of options
         // that can be changed from the
         // parameter file.
         //
         // There are two more practicalities:
         // Since we have built our right hand
         // side and solution vector as
         // deal.II Vector objects (as opposed
         // to the matrix, which is a Trilinos
         // object), we must hand the solvers
         // Trilinos Epetra vectors.  Luckily,
         // they support the concept of a
         // 'view', so we just send in a
         // pointer to our deal.II vectors. We
         // have to provide an Epetra_Map for
         // the vector that sets the parallel
         // distribution, which is just a
         // dummy object in serial. The
         // easiest way is to ask the matrix
         // for its map, and we're going to be
         // ready for matrix-vector products
         // with it.
         //
         // Secondly, the Aztec solver wants
         // us to pass a Trilinos
         // Epetra_CrsMatrix in, not the
         // deal.II wrapper class itself. So
         // we access to the actual Trilinos
         // matrix in the Trilinos wrapper
         // class by the command
         // trilinos_matrix(). Trilinos wants
         // the matrix to be non-constant, so
         // we have to manually remove the
         // constantness using a const_cast.
      case Parameters::Solver::gmres:
      {
         Epetra_Vector x(View, system_matrix.domain_partitioner(),
                         newton_update.begin());
         Epetra_Vector b(View, system_matrix.range_partitioner(),
                         right_hand_side.begin());
         
         AztecOO solver;
         solver.SetAztecOption(AZ_output,
                               (parameters.output ==
                                Parameters::Solver::quiet
                                ?
                                AZ_none
                                :
                                AZ_all));
         solver.SetAztecOption(AZ_solver, AZ_gmres);
         solver.SetRHS(&b);
         solver.SetLHS(&x);
         
         
         solver.SetAztecOption(AZ_precond,         AZ_dom_decomp);
         solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
         solver.SetAztecOption(AZ_overlap,         0);
         solver.SetAztecOption(AZ_reorder,         0);
         
         solver.SetAztecParam(AZ_drop,      parameters.ilut_drop);
         solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill);
         solver.SetAztecParam(AZ_athresh,   parameters.ilut_atol);
         solver.SetAztecParam(AZ_rthresh,   parameters.ilut_rtol);
         
         solver.SetUserMatrix(const_cast<Epetra_CrsMatrix*>
                              (&system_matrix.trilinos_matrix()));
         
         solver.Iterate(parameters.max_iterations, parameters.linear_residual);
         
         return std::pair<unsigned int, double> (solver.NumIters(),
                                                 solver.TrueResidual());
      }

      /*
      case Parameters::Solver::dgmres:
      {
         SolverControl   solver_control (parameters.max_iterations, 
                                         parameters.linear_residual);
         SolverGMRES<>   solver (solver_control);

         PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;
         preconditioner.initialize(system_matrix, fe.dofs_per_cell);

         solver.solve (system_matrix, 
                       newton_update, 
                       right_hand_side,
                       preconditioner);
      }
      */
   }
   
   Assert (false, ExcNotImplemented());
   return std::pair<unsigned int, double> (0,0);
}


// @sect4{ConservationLaw::compute_refinement_indicators}

// This function is real simple: We don't
// pretend that we know here what a good
// refinement indicator would be. Rather, we
// assume that the <code>EulerEquation</code>
// class would know about this, and so we
// simply defer to the respective function
// we've implemented there:
template <int dim>
void
ConservationLaw<dim>::
compute_refinement_indicators (Vector<double> &refinement_indicators) const
{
   EulerEquations<dim>::compute_refinement_indicators (dof_handler,
                                                       mapping,
                                                       predictor,
                                                       refinement_indicators);
}



// @sect4{ConservationLaw::refine_grid}

// Here, we use the refinement indicators
// computed before and refine the mesh. At
// the beginning, we loop over all cells and
// mark those that we think should be
// refined:
template <int dim>
void
ConservationLaw<dim>::refine_grid (const Vector<double> &refinement_indicators)
{
   typename DoFHandler<dim>::active_cell_iterator
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();
   
   for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
   {
      cell->clear_coarsen_flag();
      cell->clear_refine_flag();
      
      if ((cell->level() < parameters.shock_levels) &&
          (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val))
         cell->set_refine_flag();
      else
         if ((cell->level() > 0) &&
             (std::fabs(refinement_indicators(cell_no)) < 0.75*parameters.shock_val))
            cell->set_coarsen_flag();
   }
   
   // Then we need to transfer the
   // various solution vectors from
   // the old to the new grid while we
   // do the refinement. The
   // SolutionTransfer class is our
   // friend here; it has a fairly
   // extensive documentation,
   // including examples, so we won't
   // comment much on the following
   // code. The last three lines
   // simply re-set the sizes of some
   // other vectors to the now correct
   // size:
   std::vector<Vector<double> > transfer_in;
   std::vector<Vector<double> > transfer_out;
   
   transfer_in.push_back(old_solution);
   transfer_in.push_back(predictor);
   
   triangulation.prepare_coarsening_and_refinement();
   
   SolutionTransfer<dim> soltrans(dof_handler);
   soltrans.prepare_for_coarsening_and_refinement(transfer_in);
   
   triangulation.execute_coarsening_and_refinement ();
   
   dof_handler.clear();
   dof_handler.distribute_dofs (fe);
   
   {
      Vector<double> new_old_solution(1);
      Vector<double> new_predictor(1);
      
      transfer_out.push_back(new_old_solution);
      transfer_out.push_back(new_predictor);
      transfer_out[0].reinit(dof_handler.n_dofs());
      transfer_out[1].reinit(dof_handler.n_dofs());
   }
   
   soltrans.interpolate(transfer_in, transfer_out);
   
   old_solution.reinit (transfer_out[0].size());
   old_solution = transfer_out[0];
   
   predictor.reinit (transfer_out[1].size());
   predictor = transfer_out[1];
   
   current_solution.reinit(dof_handler.n_dofs());
   current_solution = old_solution;
   right_hand_side.reinit (dof_handler.n_dofs());
}


// @sect4{ConservationLaw::output_results}

// This function now is rather
// straightforward. All the magic, including
// transforming data from conservative
// variables to physical ones has been
// abstracted and moved into the
// EulerEquations class so that it can be
// replaced in case we want to solve some
// other hyperbolic conservation law.
//
// Note that the number of the output file is
// determined by keeping a counter in the
// form of a static variable that is set to
// zero the first time we come to this
// function and is incremented by one at the
// end of each invokation.
template <int dim>
void ConservationLaw<dim>::output_results () const
{
   typename EulerEquations<dim>::Postprocessor
   postprocessor (parameters.schlieren_plot);
   
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   
   data_out.add_data_vector (current_solution,
                             EulerEquations<dim>::component_names (),
                             DataOut<dim>::type_dof_data,
                             EulerEquations<dim>::component_interpretation ());
   
   data_out.add_data_vector (current_solution, postprocessor);
   
   data_out.build_patches (fe.degree);
   
   static unsigned int output_file_number = 0;
   std::string filename = "solution-" + Utilities::int_to_string (output_file_number, 3);
     
   if(parameters.output_format == "vtk")     
      filename += ".vtk";
   else if(parameters.output_format == "tecplot") 
      filename += ".plt";

   std::ofstream output (filename.c_str());

   if(parameters.output_format == "vtk")     
      data_out.write_vtk (output);
   else if(parameters.output_format == "tecplot") 
      data_out.write_tecplot (output);
   
   ++output_file_number;
}




// @sect4{ConservationLaw::run}

// This function contains the top-level logic
// of this program: initialization, the time
// loop, and the inner Newton iteration.
//
// At the beginning, we read the mesh file
// specified by the parameter file, setup the
// DoFHandler and various vectors, and then
// interpolate the given initial conditions
// on this mesh. We then perform a number of
// mesh refinements, based on the initial
// conditions, to obtain a mesh that is
// already well adapted to the starting
// solution. At the end of this process, we
// output the initial solution.
template <int dim>
void ConservationLaw<dim>::run ()
{
   {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      
      std::ifstream input_file(parameters.mesh_filename.c_str());
      Assert (input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));
      
      if(parameters.mesh_type == "ucd")
         grid_in.read_ucd(input_file);
      else if(parameters.mesh_type == "gmsh")
         grid_in.read_msh(input_file);
   }

   /*
   static const HyperBallBoundary<dim> boundary_description;
   triangulation.set_boundary (1, boundary_description);
   */
   
   dof_handler.clear();
   dof_handler.distribute_dofs (fe);
   
   // Size all of the fields.
   old_solution.reinit (dof_handler.n_dofs());
   current_solution.reinit (dof_handler.n_dofs());
   predictor.reinit (dof_handler.n_dofs());
   right_hand_side.reinit (dof_handler.n_dofs());
   
   setup_system();
   
   VectorTools::interpolate(dof_handler,
                            parameters.initial_conditions, old_solution);
   current_solution = old_solution;
   predictor = old_solution;
   
   if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         
         compute_refinement_indicators(refinement_indicators);
         refine_grid(refinement_indicators);
         
         setup_system();
         
         VectorTools::interpolate(dof_handler,
                                  parameters.initial_conditions, old_solution);
         current_solution = old_solution;
         predictor = old_solution;
      }
   
   output_results ();
   
   // We then enter into the main time
   // stepping loop. At the top we simply
   // output some status information so one
   // can keep track of where a computation
   // is, as well as the header for a table
   // that indicates progress of the nonlinear
   // inner iteration:
   Vector<double> newton_update (dof_handler.n_dofs());
   
   double time = 0;
   double next_output = time + parameters.output_step;
   
   predictor = old_solution;
   while (time < parameters.final_time)
   {
      std::cout << "T=" << time << std::endl
		          << "   Number of active cells:       "
		          << triangulation.n_active_cells()
		          << std::endl
		          << "   Number of degrees of freedom: "
		          << dof_handler.n_dofs()
		          << std::endl
		          << std::endl;
      
      std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                << "   _____________________________________" << std::endl;
      
      // Then comes the inner Newton
      // iteration to solve the nonlinear
      // problem in each time step. The way
      // it works is to reset matrix and
      // right hand side to zero, then
      // assemble the linear system. If the
      // norm of the right hand side is small
      // enough, then we declare that the
      // Newton iteration has
      // converged. Otherwise, we solve the
      // linear system, update the current
      // solution with the Newton increment,
      // and output convergence
      // information. At the end, we check
      // that the number of Newton iterations
      // is not beyond a limit of 10 -- if it
      // is, it appears likely that
      // iterations are diverging and further
      // iterations would do no good. If that
      // happens, we throw an exception that
      // will be caught in
      // <code>main()</code> with status
      // information being displayed before
      // the program aborts.
      //
      // Note that the way we write the
      // AssertThrow macro below is by and
      // large equivalent to writing
      // something like <code>if
      // (!(nonlin_iter @<= 10)) throw
      // ExcMessage ("No convergence in
      // nonlinear solver");</code>. The only
      // significant difference is that
      // AssertThrow also makes sure that the
      // exception being thrown carries with
      // it information about the location
      // (file name and line number) where it
      // was generated. This is not overly
      // critical here, because there is only
      // a single place where this sort of
      // exception can happen; however, it is
      // generally a very useful tool when
      // one wants to find out where an error
      // occurred.
      unsigned int nonlin_iter = 0;
      double res_norm0 = 1.0;
      
      //if(parameters.is_stationary == false)
      current_solution = predictor;
      
      Integrator<dim> integrator (dof_handler);
      setup_mesh_worker (integrator);
      while (nonlin_iter < parameters.max_nonlin_iter)
      {
         system_matrix = 0;
         
         right_hand_side = 0;
         assemble_system (integrator);
         
         const double res_norm = right_hand_side.l2_norm();
         if(nonlin_iter == 0) res_norm0 = res_norm;
         if (std::fabs(res_norm) < 1.0e-10)
         {
            std::printf("   %-16.3e (converged)\n\n", res_norm);
            break;
         }
         else
         {
            newton_update = 0;
            
            std::pair<unsigned int, double> convergence
            = solve (newton_update);
            
            current_solution += newton_update;
            
            std::printf("   %-16.3e %04d        %-5.2e\n",
                        res_norm, convergence.first, convergence.second);
         }
         
         ++nonlin_iter;
         AssertThrow (nonlin_iter <= parameters.max_nonlin_iter, 
                      ExcMessage ("No convergence in nonlinear solver"));
      }
      
      // We only get to this point if the
      // Newton iteration has converged, so
      // do various post convergence tasks
      // here:
      //
      // First, we update the time
      // and produce graphical output
      // if so desired. Then we
      // update a predictor for the
      // solution at the next time
      // step by approximating
      // $\mathbf w^{n+1}\approx
      // \mathbf w^n + \delta t
      // \frac{\partial \mathbf
      // w}{\partial t} \approx
      // \mathbf w^n + \delta t \;
      // \frac{\mathbf w^n-\mathbf
      // w^{n-1}}{\delta t} = 2
      // \mathbf w^n - \mathbf
      // w^{n-1}$ to try and make
      // adaptivity work better.  The
      // idea is to try and refine
      // ahead of a front, rather
      // than stepping into a coarse
      // set of elements and smearing
      // the old_solution.  This
      // simple time extrapolator
      // does the job. With this, we
      // then refine the mesh if so
      // desired by the user, and
      // finally continue on with the
      // next time step:
      time += parameters.time_step;

      //parameters.time_step *= 2;
      
      if (parameters.output_step < 0)
         output_results ();
      else if (time >= next_output)
      {
         output_results ();
         next_output += parameters.output_step;
      }
      
      predictor = current_solution;
      predictor.sadd (2.0, -1.0, old_solution);
      
      old_solution = current_solution;
      
      if (parameters.do_refine == true)
      {
         Vector<double> refinement_indicators (triangulation.n_active_cells());
         compute_refinement_indicators(refinement_indicators);
         
         refine_grid(refinement_indicators);
         setup_system();
         
         newton_update.reinit (dof_handler.n_dofs());
      }
   }
}

template class ConservationLaw<2>;
