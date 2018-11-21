#ifndef SOLUTION_OPERATIONS__
#define SOLUTION_OPERATIONS__

#include "general.h"

//polymorphic base class to define additional solution operations from (transfer, update, etc.)
template <class Type>
class SolutionOperation
{
public:
  //pure virtual members, must be overloaded in derived operations class
  virtual ~SolutionOperation(){};
  virtual void Apply() = 0;
private:
};


#endif
