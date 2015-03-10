#ifndef __FPARSER_H__
#define __FPARSER_H__

#include "fparser.hh"
#include <cmath>
#include <map>

extern std::map<std::string, double> constants;

class FParser: public FunctionParser
{
   public:
      FParser()
      {
         AddConstant("pi", M_PI);
      }
      void FParse (const std::string& fun)
      {
         std::map<std::string,double>::iterator it;
         for ( it=constants.begin() ; it != constants.end(); it++ )
            AddConstant ((*it).first, (*it).second);
         Parse (fun, "x,y");
      }
};

#endif
