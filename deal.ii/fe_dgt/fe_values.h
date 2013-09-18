  /** Computes shape functions, gradients, etc. at the specified
   *  @p points. We do not compute quadrature points, jacobians, normals
   *  etc. which depend on the cell mapping. This will be used only
   *  while computing the WENO limiter.
   */
  template <class DH, bool level_dof_access>
  void reinit 
     (
     const TriaIterator<DoFCellAccessor<DH,level_dof_access> > cell,
     const std::vector< Point<spacedim> >& points
     );

  /** Computes shape functions, gradients, etc. at the specified
   *  @p points. We do not compute quadrature points, jacobians, normals
   *  etc. which depend on the cell mapping. This will be used only
   *  while computing the WENO limiter.
   */
  void reinit 
     (
     const typename Triangulation<dim,spacedim>::cell_iterator &cell,
     const std::vector< Point<spacedim> >& points
     );
