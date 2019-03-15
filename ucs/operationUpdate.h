#ifndef SOLUTION_UPDATE__
#define SOLUTION_UPDATE__

//class which defines the solution update operation
//#include "solutionOperations.h"

template <class Type>
class OperationUpdate : public SolutionOperation<Type>
{
public:
  OperationUpdate(SolutionSpaceBase<Type>& space):
    space(space)
  {};
  ~OperationUpdate(){};
  void Apply(){
    bool isConverged = space.NewtonIterate();
  };
private:
  SolutionSpaceBase<Type>& space;
};

#endif
