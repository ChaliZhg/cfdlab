#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define SIGN(a) (((a)<0) ? -1:1)

const double GAMMA = 1.4;
const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};

enum Model {euler, ns};

using namespace std;

//------------------------------------------------------------------------------
// Minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a,
               const double& b,
               const double& c)
{
   double result;
   
   if (a*b > 0.0 && b*c > 0.0)
   {
      result  = min( min(fabs(a), fabs(b)), fabs(c) );
      result *= SIGN(a);
   }
   else
      result = 0.0;
   
   return result;
   
}

//------------------------------------------------------------------------------
// Reconstruct left state of right face
//------------------------------------------------------------------------------
vector<double> muscl (const vector<double>& ul,
                      const vector<double>& uc,
                      const vector<double>& ur)
{
   unsigned int n = ul.size();
   vector<double> result (n);
   double dul, duc, dur;
   
   for(unsigned int i=0; i<n; ++i)
   {
      dul = uc[i] - ul[i];
      dur = ur[i] - uc[i];
      duc = (ur[i] - ul[i])/2.0;
      result[i] = uc[i] + 0.5 * minmod (1.5*dul, duc, 1.5*dur);
   }
   
   return result;
}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
class FVProblem
{
   public:
      FVProblem ();
      void run ();
   
   private:
      void make_grid_and_dofs ();
      void initial_condition ();
      void compute_dt ();
      void con_to_prim ();
      void compute_face_derivatives ();
      void reconstruct (const unsigned int face,
                        vector<double>& left,
                        vector<double>& right) const;
      void kfvs_split_flux (const vector<double>& prim,
                            const int sign,
                            vector<double>& flux) const;
      void num_flux (const vector<double>&,
                     const vector<double>&,
                     vector<double>&) const;
      void compute_residual ();
      void update_solution (const unsigned int rk);
      void output ();
      
      double d_left, u_left, p_left;
      double d_right, u_right, p_right;
      unsigned int n_var;
      unsigned int n_cell;
      unsigned int n_face;
      double xmin;
      double xmax;
      double dx;
      double dt;
      double cfl;
      double final_time;
      Model  model;
      vector<double> xc;
      vector<double> xf;
   
      vector< vector<double> > primitive;
      vector< vector<double> > residual;
      vector< vector<double> > conserved;
      vector< vector<double> > conserved_old;
      vector<double> ux_f, Tx_f;

};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
FVProblem::FVProblem ()
{   
   model  = euler;
   n_var  = 3;
   n_cell = 100;
   xmin   = 0.0;
   xmax   = 1.0;
   cfl    = 0.8;
   final_time = 0.2;
   
   d_left  = 1.0;
   d_right = 0.125;
   
   u_left  = 0.0;
   u_right = 0.0;
   
   p_left  = 1.0;
   p_right = 0.1;
}

//------------------------------------------------------------------------------
// Allocate memory for grid and create grid
//------------------------------------------------------------------------------
void FVProblem::make_grid_and_dofs ()
{
   cout << "Making grid and allocating memory ...\n";
   
   n_face = n_cell + 1;
   dx = (xmax - xmin) / n_cell;
   xc.resize (n_cell);
   xf.resize (n_face);
   
   // Make grid
   for(unsigned int i=0; i<n_face; ++i)
      xf[i] = xmin + i * dx;
   for(unsigned int i=0; i<n_cell; ++i)
      xc[i] = 0.5 * (xf[i] + xf[i+1]);
   
   primitive.resize (n_cell, vector<double>(n_var));
   residual.resize (n_cell, vector<double>(n_var));
   conserved.resize (n_cell, vector<double>(n_var));
   conserved_old.resize (n_cell, vector<double>(n_var));
   
   if(model == ns)
   {
      ux_f.resize (n_face);
      Tx_f.resize (n_face);
   }
}

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FVProblem::initial_condition ()
{
   cout << "Setting initial conditions ...\n";
   
   // Set initial condition
   for(unsigned int i=0; i<n_cell; ++i)
   {
      if(xf[i+1] <= 0.5)
      {
         conserved[i][0] = d_left;
         conserved[i][1] = d_left * u_left;
         conserved[i][2] = p_left/(GAMMA-1.0) + 0.5 * d_left * pow(u_left,2);
      }
      else
      {
         conserved[i][0] = d_right;
         conserved[i][1] = d_right * u_right;
         conserved[i][2] = p_right/(GAMMA-1.0) + 0.5 * d_right * pow(u_right,2);
      }
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FVProblem::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double speed = fabs(primitive[i][1]) + 
                     sqrt(GAMMA * primitive[i][2] / primitive[i][0]);
      dt = min (dt, 1.0/speed);
   }
   dt *= cfl * dx;
}

// Convert conserved to primitive
void FVProblem::con_to_prim ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
      primitive[i][0] = conserved[i][0];
      primitive[i][1] = conserved[i][1]/conserved[i][0];
      primitive[i][2] = (GAMMA-1.0) * (conserved[i][2] - 
                           0.5 * pow(conserved[i][1], 2.0) / conserved[i][0]);
   }
}

//------------------------------------------------------------------------------
// Compute derivatives at faces
//------------------------------------------------------------------------------
void FVProblem::compute_face_derivatives ()
{
   ux_f[0] = 0.0;
   Tx_f[0] = 0.0;
   
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      ux_f[i] = (primitive[i][1] - primitive[i-1][1]) / dx;
      Tx_f[i] = (primitive[i][2]   / primitive[i][0] -
                 primitive[i-1][2] / primitive[i-1][0]) / dx;
   }
   
   ux_f[n_face-1] = 0.0;
   Tx_f[n_face-1] = 0.0;
}

//------------------------------------------------------------------------------
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FVProblem::reconstruct (const unsigned int face,
                             vector<double>& left,
                             vector<double>& right) const
{
   if(face==1)
   {
      left = primitive[0];
      right = muscl (primitive[2], primitive[1], primitive[0]);
   }
   else if(face==n_face-2)
   {
      left  = muscl (primitive[n_cell-3], primitive[n_cell-2], 
                     primitive[n_cell-1]);
      right = primitive[n_cell-1];
   }
   else
   {
      left  = muscl (primitive[face-2], primitive[face-1], primitive[face]);
      right = muscl (primitive[face+1], primitive[face], primitive[face-1]);
   }
}

//------------------------------------------------------------------------------
// KFVS split fluxes: sign=+1 give positive flux and
// sign=-1 gives negative flux
//------------------------------------------------------------------------------
void FVProblem::kfvs_split_flux (const vector<double>& prim,
                                 const int sign,
                                 vector<double>& flux) const
{
   double beta, s, A, B, E, fact;

   beta = 0.5 * prim[0] / prim[2];
   s    = prim[1] * sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   E    = prim[2]/(GAMMA-1.0) + 0.5 * prim[0] * pow(prim[1], 2);
   fact = prim[1] * A + B;
   
   flux[0] = prim[0] * fact;
   flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A + 
             prim[0] * prim[1] * B;
   flux[2] = prim[1] * (E + prim[2]) * A +
             (E + 0.5 * prim[2]) * B;
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FVProblem::num_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>&       flux) const
{
   vector<double> flux_pos(n_var);
   vector<double> flux_neg(n_var);
   
   kfvs_split_flux (left,  +1, flux_pos);
   kfvs_split_flux (right, -1, flux_neg);
   
   for(unsigned int i=0; i<n_var; ++i)
      flux[i] = flux_pos[i] + flux_neg[i];
   
  
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         residual[i][j] = 0.0;
   
   vector<double> flux (n_var);
   vector<double> left (n_var);
   vector<double> right(n_var);
   
   // Flux through left boundary
   num_flux (primitive[0], primitive[0], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[0][j] -= flux[j];

   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (left, right, flux);
      for(unsigned int j=0; j<n_var; ++j)
      {
         residual[i-1][j] += flux[j];
         residual[i][j]   -= flux[j];
      }
   }

   // Flux through right boundary
   num_flux (primitive[n_cell-1], primitive[n_cell-1], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[n_cell-1][j] += flux[j];
}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::update_solution (const unsigned int rk)
{
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         conserved[i][j] = arks[rk] * conserved_old[i][j] +
            brks[rk] * (conserved[i][j] - (dt/dx) * residual[i][j]);
   
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FVProblem::output ()
{
   cout << "Saving solution to sol.dat\n";
   
   ofstream fo("sol.dat");
   for(unsigned int i=0; i<n_cell; ++i)
      fo << xc[i] << " " 
         << primitive[i][0] << " " 
         << primitive[i][1] << " "
         << primitive[i][2] << endl;
   fo.close ();
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FVProblem::run ()
{
   make_grid_and_dofs ();
   initial_condition ();
   con_to_prim ();

   double time = 0.0;
   unsigned int iter = 0;
   while (time < final_time)
   {
      conserved_old = conserved;
      compute_dt ();
      if(time+dt > final_time) dt = final_time - time;
      for(unsigned int rk=0; rk<3; ++rk)
      {
         compute_residual ();
         update_solution (rk);
         con_to_prim ();
      }
      time += dt;
      ++iter;
      cout << "Iter = " << iter << " Time = " << time << endl;
   }
   
   con_to_prim ();
   output ();
}

//------------------------------------------------------------------------------
int main ()
{
   FVProblem fv_problem;
   fv_problem.run ();
   
   return 0;
}
