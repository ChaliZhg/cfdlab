#ifndef __BC_H__
#define __BC_H__

#include <iostream>
#include <string>
#include <cassert>
#include "face.h"
#include "fparser.h"
#include "primvar.h"
#include "material.h"

namespace BC
{
   enum BCType { none, slip, noslip, maxwell, farfield, inlet, outlet, pressure };
}

//------------------------------------------------------------------------------
// Boundary condition class
//------------------------------------------------------------------------------
class BoundaryCondition
{
   public:
      BoundaryCondition () {};
      BoundaryCondition (Material                 &material,
                         std::string              &bc_type,
                         std::vector<std::string> &variable,
                         std::vector<std::string> &function);
      void apply (const Vector        &vertex,
                  const Face          &face,
                  std::vector<ConVar> &state);
      void apply (const Vector &vertex,
                  const Face   &face,
                  ConVar       &state);
      void apply_slip (const Face          &face,
                       std::vector<ConVar> &state);
      void apply_slip (const Face &face,
                       ConVar     &state);
      void apply_noslip (const Vector        &vertex,
                         std::vector<ConVar> &state);
      void apply_noslip (const Vector &vertex,
                         ConVar       &state);
      void apply_maxwell (const Face          &face,
                          std::vector<ConVar> &state);
      void apply_pressure (const Vector        &vertex,
                           std::vector<ConVar> &state);
      void apply_inlet (const Vector        &vertex,
                        std::vector<ConVar> &state);
      void apply_inlet (const Vector &vertex,
                        ConVar       &state);
      void apply_outlet (std::vector<ConVar> &state);
      void apply_farfield (const Vector        &vertex,
                           std::vector<ConVar> &state);
      std::string    name;
      BC::BCType     type;
      bool           adiabatic;

   private:
      Material*      material;
      FParser        xvelocity;
      FParser        yvelocity;
      FParser        zvelocity;
      FParser        pressure;
      FParser        density;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
inline
BoundaryCondition::BoundaryCondition (Material                 &material,
                                      std::string              &bc_type,
                                      std::vector<std::string> &variable,
                                      std::vector<std::string> &function)
:
   name (bc_type),
   adiabatic (false),
   material(&material)
{
   // Set to none for safety purpose
   type = BC::none;

   // Slip bc, no state is required
   if(bc_type == "slip")
   {
      assert (variable.size() == 0);
      type = BC::slip;
      adiabatic = true;
   }
   // noslip bc: velocity is specified. If temperature is also specified
   // then it is also used. In this case, adiabatic bc is not used
   else if(bc_type == "noslip")
   {
      assert (variable.size() == 3 || variable.size() == 4);
      type = BC::noslip;
      adiabatic = true;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.FParse (function[i]);
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.FParse (function[i]);
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.FParse (function[i]);
         }
         else if(variable[i] == "density")
         {
            density.FParse (function[i]);
         }
      }
      assert (has_xvelocity && has_yvelocity && has_zvelocity);
   }
   else if(bc_type == "maxwell")
   {
      assert (variable.size() == 4);
      type = BC::maxwell;
      adiabatic = false;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      bool has_density = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.FParse (function[i]);
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.FParse (function[i]);
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.FParse (function[i]);
         }
         else if(variable[i] == "density")
         {
            has_density = true;
            density.FParse (function[i]);
         }
      }
      assert (has_xvelocity && has_yvelocity && has_zvelocity && has_density);
   }
   // In this case only pressure is specified
   else if(bc_type == "pressure")
   {
      assert (variable.size() == 1);
      type = BC::pressure;
      assert (variable[0] == "pressure");
      pressure.FParse (function[0]);
   }
   // All values are specified
   else if(bc_type == "inlet" || bc_type == "farfield")
   {
      assert (variable.size() == 5);
      if(bc_type == "inlet")
         type = BC::inlet;
      else
         type = BC::farfield;
      bool has_density = false;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      bool has_pressure    = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "density")
         {
            has_density = true;
            density.FParse (function[i]);
         }
         else if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.FParse (function[i]);
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.FParse (function[i]);
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.FParse (function[i]);
         }
         else if(variable[i] == "pressure")
         {
            has_pressure = true;
            pressure.FParse (function[i]);
         }
      }
      assert (has_density && has_xvelocity && has_yvelocity && has_zvelocity &&
              has_pressure);
   }
   // At outflow nothing is specified
   else if(bc_type == "outlet")
   {
      assert (variable.size() == 0);
      type = BC::outlet;
   }
   else
   {
      std::cout << "BoundaryCondition: Unknown boundary condition " << bc_type << std::endl;
      abort();
   }

   if(type == BC::none)
   {
      std::cout << "BoundaryCondition: unknown bc for " << bc_type << std::endl;
      abort();
   }
}

//------------------------------------------------------------------------------
// Normal velocity is zero
// u[1] = u[0] - 2(u[0] * n) n
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_slip(const Face          &face,
                                   std::vector<ConVar> &state)
{
   Vector unit_normal = face.normal / face.measure;
   state[1].momentum -= unit_normal * (state[1].momentum * unit_normal) * 2.0;
}

//------------------------------------------------------------------------------
// At inlet all values are specified
// Both states are set to inlet values
// e.g., supersonic inlet
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_inlet (const Vector         &vertex,
                                     std::vector<ConVar> &state)
{
   double point[2]  = {vertex.x, vertex.y};
   double rho = density.Eval(point);
   double pre = pressure.Eval(point);

   Vector velocity;
   velocity.x = xvelocity.Eval(point);
   velocity.y = yvelocity.Eval(point);
   
   state[0].density  = rho;
   state[0].momentum = velocity * rho;
   state[0].energy   = pre/(material->gamma-1) + 0.5 * rho * velocity.square();
   
   state[1] = state[0];
}

//------------------------------------------------------------------------------
// At outlet all all values are from inside values
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_outlet (std::vector<ConVar> &state)
{
   state[1] = state[0];
}

//------------------------------------------------------------------------------
// At farfield all values are specified
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_farfield (const Vector        &vertex,
                                        std::vector<ConVar> &state)
{
   double point[2]     = {vertex.x, vertex.y};
   
   double rho = density.Eval(point);
   double pre = pressure.Eval(point);
   
   Vector velocity;
   velocity.x = xvelocity.Eval(point);
   velocity.y = yvelocity.Eval(point);
   
   state[1].density = rho;
   state[1].momentum = velocity * rho;
   state[1].energy = pre/(material->gamma-1) + 0.5 * rho * velocity.square();
}

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector        &vertex,
                              const Face          &face,
                              std::vector<ConVar> &state)
{
   switch(type)
   {
      case BC::slip:
         apply_slip (face, state);
         break;

      case BC::inlet:
         apply_inlet (vertex, state);
         break;

      case BC::outlet:
         apply_outlet (state);
         break;

      case BC::farfield:
         apply_farfield (vertex, state);
         break;

      default:
         std::cout << "BoundaryCondition::apply" << std::endl;
         std::cout << "   Unknown boundary condition: " << name << std::endl;
         abort ();
   }
}

#endif
