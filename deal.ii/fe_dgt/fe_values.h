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
