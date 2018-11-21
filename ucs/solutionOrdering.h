#ifndef SOLUTION_ORDERING_H__
#define SOLUTION_ORDERING_H__

#include <vector>
#include "solutionOperations.h"
#include "solutionSpaceBase.h"
#include "operationTransfer.h"
#include "operationUpdate.h"

template <class Type>
class SolutionOrdering
{
public:
  template <class Type2>
  SolutionOrdering<Type>& operator= (const SolutionOrdering<Type2>& orderingToCopy); 
  ~SolutionOrdering();
  UInt size();
  Int Finalize(std::vector<SolutionSpaceBase<Type>*>& solSpaces);
  void Print();
  void Iterate();
  Int Insert(const std::string lineCommand);
  Int Read(std::string casename, std::string pathname);
  std::vector<std::string> commands;
private:
  Int ReadOrderingSegment(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd);
  std::vector<SolutionOperation<Type>*> opList;
};

//include implementations
#include "solutionOrdering.tcc"

#endif
