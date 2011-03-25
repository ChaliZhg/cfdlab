#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

/*
 *  integrator.h
 *  
 *
 *  Created by Praveen Chandrashekar on 25/02/11.
 *  Copyright 2011 TIFR-CAM. All rights reserved.
 *
 */

#include <dofs/dof_handler.h>

#include <lac/vector.h>
#include <lac/trilinos_sparse_matrix.h>

#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_info.h>
#include <numerics/mesh_worker_assembler.h>
#include <numerics/mesh_worker_loop.h> 


//------------------------------------------------------------------------------
// Class for integrating rhs using MeshWorker
//------------------------------------------------------------------------------
template <int dim>
class Integrator
{
public:
   Integrator (const dealii::DoFHandler<dim>& dof_handler)
      : 
      dof_info (dof_handler) {};
   
   dealii::MeshWorker::IntegrationInfoBox<dim> info_box;
   dealii::MeshWorker::DoFInfo<dim> dof_info;
   dealii::MeshWorker::Assembler::SystemSimple<dealii::TrilinosWrappers::SparseMatrix, dealii::Vector<double> > assembler;

};

#endif
