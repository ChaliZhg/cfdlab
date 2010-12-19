#ifndef __GRID_H__
#define __GRID_H__

#include <vector>

#define NEW    0
#define OLD    1

// Stochastic sample
template <int dim>
class Sample
{
   public:
      Sample (const unsigned int n_var,
              const unsigned int n_cell,
              const unsigned int counter);
      ~Sample () {};
      void clear ();
      void read ();

      double x[dim];
      unsigned int n_var, n_cell;
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
               const unsigned int n_moment);
      ~Element () {};
   
      int order;                // order = linear (1) or quadratic (2)
      unsigned int n_moment;
      unsigned int n_dof;       // no. of dof for this element
      // pointer to sample for each dof
      std::vector<typename Sample<dim>::Sample*> dof; 
      int status;               // NEW or OLD.
      Element* parent;          // parent element
      double* moment;           // element contribution to moment
      double* adj_cor;          // element contribution to adjoint correction
      bool active;

};

// Stochastic grid
template <int dim>
class Grid
{
   public:
      Grid () { n_element = 0; };
      ~Grid () {};

      unsigned int n_element;
      std::vector<typename Element<dim>::Element> element;
};

#endif
