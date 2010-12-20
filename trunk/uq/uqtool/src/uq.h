#ifndef __UQ_H__
#define __UQ_H__

#include <vector>
#include "grid.h"
#include "pdf.h"

#define STOCHASTIC   1
#define COMBINED     2

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
      void flag_elements ();
      void refine_grid ();
      void output (const unsigned int iter) const;

      PDFData<dim> pdf_data;
   
      unsigned int n_moment;
      double* moment;
      double* adj_cor;
      double* RE;

      unsigned int n_sample;
      std::vector<typename Sample<dim>::Sample> sample;
      unsigned int max_sample;

      Grid<dim> grid;
   
      unsigned int n_var; // no. of variables per cell
      unsigned int n_cell;// no. of cells
   
      int refine_type;    // STOCHASTIC or COMBINED
      double* mesh_error; // Indicator for physical mesh
};

#endif
