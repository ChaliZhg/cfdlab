#include <iostream>
#include <fstream>

using namespace std;

bool check_for_stop_file ()
{
   ifstream ifile("STOP");
   if(ifile)
      cout << "Found STOP file !!!\n";
   return ifile;
}
