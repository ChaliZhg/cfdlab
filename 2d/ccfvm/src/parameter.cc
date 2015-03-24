#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include "parameter.h"

using namespace std;
extern map<string,double> constants;

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   cout << "Reading input file " << file << endl;
   Reader fin(file);

   read_grid (fin);
   read_numeric (fin);
   read_material (fin);
   read_constants (fin);
   read_initial_condition (fin);
   read_boundary (fin);
   read_output (fin);
   check ();
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid (Reader &fin)
{
   cout << "  Reading grid section\n";

   string input;

   fin.begin_section ("grid");

   fin.entry ("type");
   fin >> input;
   if(input=="gmsh")
      grid_type = gmsh;
   else if(input=="bamg")
      grid_type = bamg;
   else if(input=="delaundo")
      grid_type = delaundo;
   else
   {
      cout << "   Unknown file type " << input << endl;
      exit (0);
   }

   fin.entry ("file");
   fin >> grid_file;

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric (Reader &fin)
{
   cout << "  Reading numeric section\n";

   string input;

   fin.begin_section ("numeric");

   fin.entry ("time_mode");
   fin >> time_mode;
   assert (time_mode == "steady" || time_mode == "unsteady");

   fin.entry ("time_scheme");
   fin >> time_scheme;
   assert (time_scheme == "rk1"    || 
           time_scheme == "ssprk3" ||
           time_scheme == "rk3"    ||
           time_scheme == "rk4"    ||
           time_scheme == "lusgs");
   if(time_scheme=="rk1")    n_rks = 1;
   if(time_scheme=="ssprk3") n_rks = 3;
   if(time_scheme=="rk3")    n_rks = 3;
   if(time_scheme=="rk4")    n_rks = 4;
   if(time_scheme=="lusgs")  n_rks = 1; 

   fin.entry ("time_step");
   fin >> time_step;
   assert (time_step >= 0.0);

   fin.entry ("cfl");
   fin >> cfl;
   assert (cfl >= 0.0);

   fin.entry ("max_iter");
   fin >> max_iter;
   assert (max_iter > 0);

   fin.entry ("final_time");
   fin >> final_time;
   assert (final_time >= 0.0);

   fin.entry ("min_residue");
   fin >> min_residue;
   assert (min_residue > 0.0);

   fin.entry ("reconstruct");
   fin >> input;
   if(input == "first")
      reconstruct_scheme = Parameter::first;
   else if(input == "second")
      reconstruct_scheme = Parameter::second;
   else if(input == "bj")
      reconstruct_scheme = Parameter::bj;
   else
   {
      cout << "read_numeric: unknown reconstruction scheme " << input << endl;
      cout << "              first, second, bj\n";
      exit (0);
   }

   fin.entry ("ang_mom");
   fin >> input;
   if(input == "yes")
      ang_mom = true;
   else if(input == "no")
      ang_mom = false;
   else
   {
      cout << "read_numeric: unknown option for ang_mom " << input << endl;
      exit (0);
   }
   
   fin.entry ("smooth_res");
   fin >> input;
   if(input == "yes")
      smooth_res = true;
   else if(input == "no")
      smooth_res = false;
   else
   {
      cout << "read_numeric: unknown option for smooth_res " << input << endl;
      exit (0);
   }

   fin.entry ("ducros");
   fin >> input;
   if(input == "yes")
      ducros = true;
   else if(input == "no")
      ducros = false;
   else
   {
      cout << "read_numeric: unknown option for ducros " << input << endl;
      exit (0);
   }

   fin.end_section ();

   // Some parameter checks
   if(time_scheme == "lusgs")
   {
      assert (time_mode != "unsteady");
      assert (time_step == 0.0);
      assert (cfl > 0.0);
   }

   if(time_step == 0.0 && cfl == 0.0)
   {
      cout << "Either time_step or cfl must be non-zero\n";
      exit (0);
   }
   if(time_step > 0.0 && cfl > 0.0)
   {
      cout << "You have specified both time_step and cfl\n";
      cout << "Specify only one of them, other being zero\n";
      exit (0);
   }

   // For steady flow, no need for final time
   if(time_mode == "steady")
      final_time = 1.0e20;
   else
   {
      assert (final_time > 0.0);
      // For unsteady flow, final time must be specified, so dont put any limit
      // on max_iter
      max_iter = 99999999;
      // For unsteady flow, dont use residual smoothing
      assert(!smooth_res);
   }

}

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material (Reader &fin)
{
   cout << "  Reading material section\n";

   string input;

   fin.begin_section ("material");

   fin.entry ("gamma");
   fin >> material.gamma;
   assert (material.gamma > 1.0);

   fin.entry ("gas_const");
   fin >> material.gas_const;
   assert (material.gas_const > 0.0);

   fin.begin_section ("viscosity");
   fin.entry ("model");
   fin >> input;
   if(input == "constant")
   {
      material.mu_model = Material::mu_constant;

      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      assert (material.mu_ref >= 0.0);
   }
   else if(input == "sutherland")
   {
      material.mu_model = Material::mu_sutherland;
      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      fin.entry ("T_ref");
      fin >> material.T_ref;
      fin.entry ("T_0");
      fin >> material.T_0;
      assert (material.mu_ref > 0.0);
      assert (material.T_ref > 0.0);
      assert (material.T_0 > 0.0);
   }
   else if(input == "power")
   {
      material.mu_model = Material::mu_power;
      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      fin.entry ("T_ref");
      fin >> material.T_ref;
      fin.entry ("omega");
      fin >> material.omega;
      assert (material.mu_ref > 0.0);
      assert (material.T_ref > 0.0);
      assert (material.omega > 0.0);
   }
   else
   {
      cout << "read_material: unknown viscosity type " << input << endl;
      exit (0);
   }
   fin.end_section ();

   fin.entry ("prandtl");
   fin >> material.prandtl;
   assert (material.prandtl > 0.0);

   fin.entry ("model");
   fin >> input;
   if(input == "euler")
      material.model = Material::euler;
   else if(input == "ns")
      material.model = Material::ns;
   else
   {
      cout << "read_material: unknown flow model " << input << endl;
      exit (0);
   }

   fin.entry ("flux");
   fin >> input;
   if(input == "kep")
      material.flux_scheme = Material::kep;
   else if(input == "lxf")
      material.flux_scheme = Material::lxf;
   else if(input == "roe")
      material.flux_scheme = Material::roe;
   else if(input == "hllc")
      material.flux_scheme = Material::hllc;
   else if(input == "kfvs")
      material.flux_scheme = Material::kfvs;
   else if(input == "kepes")
      material.flux_scheme = Material::kepes;
   else if(input == "kepes_roe")
      material.flux_scheme = Material::kepes_roe;
   else if(input == "kepes_roe2")
      material.flux_scheme = Material::kepes_roe2;
   else if(input == "kepes_rus")
      material.flux_scheme = Material::kepes_rus;
   else if(input == "kepes_hyb")
      material.flux_scheme = Material::kepes_hyb;
   else
   {
      cout << "read_material:: unknown flux scheme: " << input << endl;
      exit (0);
   }

   fin.end_section ();

   material.initialize ();
}

//------------------------------------------------------------------------------
// Read constants section
//------------------------------------------------------------------------------
void Parameter::read_constants (Reader &fin)
{
   cout << "  Reading constants section\n";

   string input;
   double value;

   fin.begin_section ("constants");

   while (!fin.eos())
   {
      fin >> input;
      fin >> value;
      cout << setw(16) << input << setw(16) << value << endl;
      constants.insert ( pair<string,double>(input, value) );
   }

}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition (Reader &fin)
{
   cout << "  Reading initial condition section\n";

   string input;

   fin.begin_section ("initial_condition");

   fin >> input;

   if(input == "library")
   {
      fin >> input;
      initial_condition.add (input);
   }
   else
   {
      assert (input == "density");
      fin.getline (input);
      initial_condition.add ("density", input);

      fin.entry ("xvelocity");
      fin.getline (input);
      initial_condition.add ("xvelocity", input);

      fin.entry ("yvelocity");
      fin.getline (input);
      initial_condition.add ("yvelocity", input);

      fin.entry ("zvelocity");
      fin.getline (input);
      initial_condition.add ("zvelocity", input);

      fin.entry ("pressure");
      fin.getline (input);
      initial_condition.add ("pressure", input);
   }

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read boundary conditions
//------------------------------------------------------------------------------
void Parameter::read_boundary (Reader &fin)
{
   cout << "  Reading boundary section\n";

   string input;

   fin.begin_section ("boundary");

   while (!fin.eos())
   {
      vector<int> f_type;
      while(!fin.bos())
      {
         int tmp;
         fin >> tmp;
         f_type.push_back (tmp);
      }

      fin.entry ("type");
      string bc_type;
      fin >> bc_type;

      vector<string> variable, function;
      while(!fin.eos())
      {
         fin >> input;
         variable.push_back (input);
         fin.getline(input);
         function.push_back (input);
      }
      BoundaryCondition bc(material, bc_type, variable, function);
      for(unsigned int i=0; i<f_type.size(); ++i)
         boundary_condition.insert (pair<int,BoundaryCondition>(f_type[i], bc));
   }
}

//------------------------------------------------------------------------------
// Read output section
//------------------------------------------------------------------------------
void Parameter::read_output (Reader &fin)
{
   cout << "  Reading output section\n";

   string input;

   fin.begin_section ("output");

   fin.entry ("format");
   fin >> write_format;
   assert (write_format == "vtk" || write_format == "tec");

   fin.entry ("frequency");
   fin >> write_frequency;
   assert (write_frequency > 0);

   fin.begin_section ("variables");
   while (!fin.eos())
   {
      fin >> input;
      assert (input=="mach" || input=="density" || input=="vorticity");
      write_variables.push_back (input);
   }

   fin.entry ("restart");
   fin >> input;
   if(input=="false")
      write_restart = false;
   else if(input=="true")
      write_restart = true;
   else
   {
      cout << "   Unknown input: " << input << endl;
      exit (0);
   }

   has_global = false;

   fin.entry ("global_KE");
   fin >> input;
   if(input=="false")
      global_KE = false;
   else if(input=="true")
   {
      global_KE = true;
      has_global= true;
   }
   else
   {
      cout << "   Unknown input: " << input << endl;
      exit (0);
   }

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Final checking of parameters
//------------------------------------------------------------------------------
void Parameter::check ()
{
   if( (material.flux_scheme == Material::kep ||
        material.flux_scheme == Material::kepes) &&
       reconstruct_scheme != Parameter::first)
   {
      reconstruct_scheme = Parameter::first;
      cout << "   kep flux is chosen; hence setting reconstruction to first order\n";
   }

}
