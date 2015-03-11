#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include <dlfcn.h>
#include "vec.h"
#include "material.h"
#include "fparser.h"

//------------------------------------------------------------------------------
// Class to store initial condition functions
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      InitialCondition ()
         : has_lib (false) {};
      void    add (std::string, std::string);
      void    add (std::string);
      PrimVar value (const Vector& p);

   private:
      bool has_lib;
      FParser density;
      FParser xvelocity;
      FParser yvelocity;
      FParser zvelocity;
      FParser pressure;

      void INITIAL_CONDITION(double,double,double&,double&,double&,double&);
      void (*initial_condition)(double,double,double&,double&,double&,double&);
};

//------------------------------------------------------------------------------
// Load function from external shared library
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string lib_file)
{
   void *lib_handle = dlopen(lib_file.c_str(), RTLD_LAZY);
   if(lib_handle == NULL)
   {
      std::cout << dlerror() << std::endl;
      abort();
   }

   initial_condition = 
      (void (*)(double,double,double&,double&,double&,double&))
      dlsym (lib_handle, "INITIAL_CONDITION");
   has_lib = true;
   
}

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string variable, std::string fun)
{
   if(variable == "density")
      density.FParse (fun);
   else if(variable == "xvelocity")
      xvelocity.FParse (fun);
   else if(variable == "yvelocity")
      yvelocity.FParse (fun);
   else if(variable == "zvelocity")
      zvelocity.FParse (fun);
   else if(variable == "pressure")
      pressure.FParse (fun);
   else
   {
      std::cout << "InitialCondition::add: Unknown variable " << variable << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
PrimVar InitialCondition::value (const Vector& p)
{
   PrimVar result;

   if (has_lib)
   {
      (*initial_condition)(p.x, p.y, 
                           result.density,
                           result.velocity.x,
                           result.velocity.y,
                           result.pressure);
   }
   else
   {
      double vals[2] = {p.x, p.y};

      result.density    = density.Eval (vals);
      result.velocity.x = xvelocity.Eval (vals);
      result.velocity.y = yvelocity.Eval (vals);
      result.pressure   = pressure.Eval (vals);
   }

   return result;
}

#endif
