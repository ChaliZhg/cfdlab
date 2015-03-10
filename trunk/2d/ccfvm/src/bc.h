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
      void apply (const Vector         &vertex,
                  const Face           &face,
                  std::vector<PrimVar> &state);
      void apply (const Vector &vertex,
                  const Face   &face,
                  PrimVar      &state);
      void apply_slip (const Face           &face,
                       std::vector<PrimVar> &state);
      void apply_slip (const Face &face,
                       PrimVar    &state);
      void apply_noslip (const Vector         &vertex,
                         std::vector<PrimVar> &state);
      void apply_noslip (const Vector &vertex,
                         PrimVar      &state);
      void apply_maxwell (const Face           &face,
                          std::vector<PrimVar> &state);
      void apply_pressure (const Vector         &vertex,
                           std::vector<PrimVar> &state);
      void apply_inlet (const Vector         &vertex,
                        std::vector<PrimVar> &state);
      void apply_inlet (const Vector &vertex,
                        PrimVar      &state);
      void apply_outlet (std::vector<PrimVar> &state);
      void apply_farfield (const Vector         &vertex,
                           std::vector<PrimVar> &state);
      std::string    name;
      BC::BCType     type;
      bool           adiabatic;

   private:
      Material*      material;
      FParser        xvelocity;
      FParser        yvelocity;
      FParser        zvelocity;
      FParser        pressure;
      FParser        temperature;
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
         else if(variable[i] == "temperature")
         {
            temperature.FParse (function[i]);
            adiabatic = false;
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
      bool has_temperature = false;
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
         else if(variable[i] == "temperature")
         {
            has_temperature = true;
            temperature.FParse (function[i]);
         }
      }
      assert (has_xvelocity && has_yvelocity && has_zvelocity && has_temperature);
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
      bool has_temperature = false;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      bool has_pressure    = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "temperature")
         {
            has_temperature = true;
            temperature.FParse (function[i]);
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
      assert (has_temperature && has_xvelocity && has_yvelocity && has_zvelocity &&
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
// u[1] = u[0] - (u[0] * n) n
// Then u[1] * n = u[0] * n = 0
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_slip(const Face           &face,
                                   std::vector<PrimVar> &state)
{
   Vector unit_normal = face.normal / face.measure;
   state[0].velocity -= unit_normal * (state[0].velocity * unit_normal);

   state[1] = state[0];
}

inline
void BoundaryCondition::apply_slip(const Face &face,
                                   PrimVar    &state)
{
   Vector unit_normal = face.normal / face.measure;
   state.velocity -= unit_normal * (state.velocity * unit_normal);
}

//------------------------------------------------------------------------------
// Velocity is specified
// If temperature is specified, then pressure is computed using the temperature
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Vector         &vertex,
                                     std::vector<PrimVar> &state)
{
   double point[2]  = {vertex.x, vertex.y};
   state[0].velocity.x = xvelocity.Eval(point);
   state[0].velocity.y = yvelocity.Eval(point);
   state[0].velocity.z = zvelocity.Eval(point);

   state[1].velocity = state[0].velocity;
   state[1].pressure = state[0].pressure;

   if(adiabatic)
      state[1].temperature  = state[0].temperature;
   else
   {
      double T = temperature.Eval(point);
      state[0].temperature = T;
      state[1].temperature = T;
   }
}

//------------------------------------------------------------------------------
// Velocity is modified, others kept same.
// Used to enforce noslip after update.
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Vector &vertex,
                                     PrimVar      &state)
{
   double point[2]  = {vertex.x, vertex.y};
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);
   if(!adiabatic)
   {
      state.temperature = temperature.Eval(point);
   }
}

//------------------------------------------------------------------------------
// state[1] = wall state
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_maxwell(const Face           &face,
                                      std::vector<PrimVar> &state)
{
   // wall values
   double point[2]  = {face.centroid.x, face.centroid.y};
   state[1].velocity.x = xvelocity.Eval(point);
   state[1].velocity.y = yvelocity.Eval(point);
   state[1].velocity.z = zvelocity.Eval(point);
   state[1].temperature = temperature.Eval(point);

   // normal velocity should be zero
   // state[1] must already have zero normal velocity
   Vector unit_normal = face.normal / face.measure;
   state[0].velocity -= unit_normal * (state[0].velocity * unit_normal);

   double density = material->Density(state[0]);

   double taunn= 0.0;
   double rhow = density * 
                 sqrt(state[0].temperature / state[1].temperature) * 
                 (1.0 - 0.5 * taunn/state[0].pressure);
   state[1].pressure = rhow * material->gas_const * state[1].temperature;
}

//------------------------------------------------------------------------------
// Reset pressure value. Other states remain same
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_pressure (const Vector         &vertex,
                                        std::vector<PrimVar> &state)
{
   double point[2]  = {vertex.x, vertex.y};

   state[1].pressure  = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At inlet all values are specified
// Both states are set to inlet values
// e.g., supersonic inlet
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_inlet (const Vector         &vertex,
                                     std::vector<PrimVar> &state)
{
   double point[2]  = {vertex.x, vertex.y};
   state[0].temperature= temperature.Eval(point);
   state[0].velocity.x = xvelocity.Eval(point);
   state[0].velocity.y = yvelocity.Eval(point);
   state[0].velocity.z = zvelocity.Eval(point);
   state[0].pressure   = pressure.Eval(point);

   state[1] = state[0];
}

inline
void BoundaryCondition::apply_inlet (const Vector &vertex,
                                     PrimVar      &state)
{
   double point[2]  = {vertex.x, vertex.y};
   state.temperature= temperature.Eval(point);
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);
   state.pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At outlet all all values are from inside values
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_outlet (std::vector<PrimVar> &state)
{
   state[1] = state[0];
}

//------------------------------------------------------------------------------
// At farfield all values are specified
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_farfield (const Vector         &vertex,
                                        std::vector<PrimVar> &state)
{
   double point[2]     = {vertex.x, vertex.y};
   state[1].temperature= temperature.Eval(point);
   state[1].velocity.x = xvelocity.Eval(point);
   state[1].velocity.y = yvelocity.Eval(point);
   state[1].velocity.z = zvelocity.Eval(point);
   state[1].pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector         &vertex,
                              const Face           &face,
                              std::vector<PrimVar> &state)
{
   switch(type)
   {
      case BC::slip:
         apply_slip (face, state);
         break;

      case BC::noslip:
         apply_noslip (vertex, state);
         break;

      case BC::maxwell:
         apply_maxwell (face, state);
         break;

      case BC::pressure:
         apply_pressure (vertex, state);
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

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector  &vertex,
                              const Face    &face,
                              PrimVar       &state)
{
   switch(type)
   {
      case BC::slip:
      case BC::maxwell:
         apply_slip (face, state);
         break;

      case BC::noslip:
         apply_noslip (vertex, state);
         break;

      case BC::pressure:
         // Nothing to be done
         break;

      case BC::inlet:
         apply_inlet (vertex, state);
         break;

      case BC::outlet:
         // Nothing to be done
         break;

      case BC::farfield:
         //apply_farfield (vertex, state);
         break;

      default:
         std::cout << "BoundaryCondition::apply" << std::endl;
         std::cout << "   Unknown boundary condition: " << name << std::endl;
         abort ();
   }
}

#endif
