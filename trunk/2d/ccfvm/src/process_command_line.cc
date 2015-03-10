#include <iostream>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include "constants.h"

extern Dimension dim;
extern bool debug;
extern bool restart;
extern bool preprocess;
extern bool bounds;
extern bool convert_to_vtk;
extern bool convert_to_tec;

using namespace std;

void show_options ();

//------------------------------------------------------------------------------
// Get command line flags and input file
//------------------------------------------------------------------------------
void process_command_line (int   argc,
                           char* argv[],
                           int&  ifile)
{ 
   if(argc < 3)
      show_options ();

   // By default, all are off
   debug      = false;
   restart    = false;
   preprocess = false;
   bounds     = false;
   convert_to_vtk = false;
   convert_to_tec = false;

   // Default is 2d
   dim        = two;

   int i = 1;
   bool found_input_file = false;

   while (i < argc)
   {
      if(strcmp(argv[i],"-axi")==0)
      {
         dim = axi;
      }
      else if(strcmp(argv[i],"-d")==0)
      {
         debug = true;
      }
      else if(strcmp(argv[i],"-r")==0)
      {
         restart = true;
      }
      else if(strcmp(argv[i],"-p")==0)
      {
         preprocess = true;
      }
      else if(strcmp(argv[i],"-b")==0)
      {
         bounds = true;
      }
      else if(strcmp(argv[i],"-vtk")==0)
      {
         convert_to_vtk = true;
      }
      else if(strcmp(argv[i],"-tec")==0)
      {
         convert_to_tec = true;
      }
      else if(strcmp(argv[i],"-i")==0)
      {
         assert (i+1 < argc); // check that there is another argument
         ifile = i+1;
         ++i;
         found_input_file = true;
      }
      else
      {
         cout << "Unknown command line flag: " << argv[i] << endl;
         show_options ();
      }

      ++i;
   }

   if(convert_to_vtk || convert_to_tec)
   {
      restart = true;
      preprocess = true;
      cout << "Conversion requested, enable reading of restart file\n";
   }

   if(!found_input_file)
      show_options ();
}

//------------------------------------------------------------------------------
// Print command line options available
//------------------------------------------------------------------------------
void show_options ()
{
   cout << "Valid flags are:\n";
   cout << "   -i filename   Specify input file name (required)\n";
   cout << "   -axi          Axisymmetric flow (optional)\n";
   cout << "   -d            Enable debug mode (optional)\n";
   cout << "   -r            Read restart file for initial condition (optional)\n";
   cout << "   -p            Do everything but do not solve (optional)\n";
   cout << "   -b            Compute min/max range of solution (optional)\n";
   cout << "   -tec          Read restart file and save in tecplot format (optional)\n";
   cout << "   -vtk          Read restart file and save in vtk format (optional)\n";
   exit (0);
}
