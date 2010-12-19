#ifndef __UQ_H__
#define __UQ_H__

#include <vector>
#include "grid.h"
#include "pdf.h"

// Main problem class
template <int dim>
class UQProblem
{
   public:
      UQProblem ();
      ~UQProblem ();
      void run ();

   private:
      void read_options ();
      void make_grid ();
      void run_simulations ();
      void compute_moments ();
      void refine_grid ();

      PDFData<dim> pdf_data;
   
      unsigned int n_moment;
      double* moment;
      double* adj_cor;

      unsigned int n_sample;
      std::vector<typename Sample<dim>::Sample> sample;
      unsigned int max_sample;

      Grid<dim> grid;
   
      unsigned int n_var; // no. of variables per cell
      unsigned int n_cell;// no. of cells
};

#endif
