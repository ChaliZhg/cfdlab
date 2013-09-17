template <int dim, int spacedim>
void FEValues<dim,spacedim>::reinit
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const std::vector< Point<spacedim> >& points)
{
   // no FE in this cell, so no assertion
   // necessary here
   this->maybe_invalidate_previous_present_cell (cell);
   this->check_cell_similarity(cell);
   
   // set new cell. auto_ptr will take
   // care that old object gets
   // destroyed and also that this
   // object gets destroyed in the
   // destruction of this class
   this->present_cell.reset
   (new typename FEValuesBase<dim,spacedim>::TriaCellIterator (cell));
   // this was the part of the work
   // that is dependent on the actual
   // data type of the iterator. now
   // pass on to the function doing
   // the real work.
   // We dont need to compute the quadrature points or mapping.
   // We directly compute the shape functions, gradients, etc.
   // at the specified points.
   
   Assert (this->quadrature_points.size() == points.size(),
           ExcDimensionMismatch(this->quadrature_points.size(), points.size()));
   for(unsigned int q=0; q<points.size(); ++q)
      this->quadrature_points[q] = points[q];
   
   this->get_fe().fill_fe_values(this->get_mapping(),
                                 *this->present_cell,
                                 quadrature,
                                 *this->mapping_data,
                                 *this->fe_data,
                                 *this,
                                 this->cell_similarity);
   
   this->fe_data->clear_first_cell ();
   this->mapping_data->clear_first_cell ();
}
