#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include <valarray>

#define NEW    0
#define OLD    1

// Stochastic sample
template <int dim>
class Sample
{
   public:
      Sample (const unsigned int n_var,
              const unsigned int n_cell,
              const unsigned int n_moment,
              const unsigned int counter);
      ~Sample () {};
      void clear ();
      void read ();

      std::valarray<double> x;
      std::valarray<double> J;
      unsigned int n_var, n_cell;
      unsigned int idx;
      double* primal;
      double* adjoint;
      char    directory[64]; // Dir in which simulation is run
      int     status;        // NEW or OLD: if simulation has been run
      bool    load;          // If primal/adjoint loaded into memory
};

// Stochastic element
template <int dim>
class Element
{
   public:
      Element (const unsigned int order,
               const unsigned int n_moment,
               const unsigned int n_cell,
               const unsigned int counter);
      ~Element () {};
      void read_dof ();
      void clear_dof ();
      void save_mesh_error ();
      void load_mesh_error ();
      std::vector<unsigned int> largest_face (const std::valarray<double>&,
                                              const std::valarray<double>&,
                                              const std::valarray<double>&);
   
      int order;             // order = linear (1) or quadratic (2)
      unsigned int n_moment; // no. of moments
      unsigned int n_cell;   // no. of physical cells
      unsigned int n_dof;    // no. of dof for this element
      // pointer to sample for each dof
      std::vector<typename Sample<dim>::Sample*> dof; 
      std::vector<unsigned int> idof;
      int status;            // NEW or OLD.
      std::valarray<double> moment;  // element contribution to moment
      std::valarray<double> adj_cor; // element contribution to adj correction
      std::valarray<double> RE;      // element contribution to remaining error
      bool active;
      bool refine_flag;
      std::valarray<double> mesh_error;
      char directory[64];

};

// Stochastic grid
template <int dim>
class Grid
{
   public:
      Grid () {};
      ~Grid () {};
      void reinit_dof (std::vector<typename Sample<dim>::Sample>& sample);
      int n_active_elements () const;

      std::vector<typename Element<dim>::Element> element;
};

#endif
