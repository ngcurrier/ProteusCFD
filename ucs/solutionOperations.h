#ifndef SOLUTION_OPERATIONS__
#define SOLUTION_OPERATIONS__

#include "general.h"
#include "transfer.h"
#include "solutionSpaceBase.h"
#include <vector>

template <class Type> class SolutionField;

template <class Type>
class SolutionOperation
{
public:
  //pure virtual members, must be overloaded in derived operations class
  virtual ~SolutionOperation(){};
  virtual void Apply() = 0;
private:
};

template <class Type>
class OperationTransfer : public SolutionOperation<Type>
{
public:
  OperationTransfer(SolutionSpaceBase<Type>& source, SolutionSpaceBase<Type>& dest, 
		    ScalarField<Type>& fieldSource, ScalarField<Type>& fieldDest);
  OperationTransfer(SolutionSpaceBase<Type>& source, SolutionSpaceBase<Type>& dest, 
		    SolutionField<Type>& fieldSource, SolutionField<Type>& fieldDest);
  ~OperationTransfer();
  void Apply();
private:
  Transfer<Type>* transfer;
  Int bcSrc;
  Int bcDest;
};

template <class Type>
class OperationUpdate : public SolutionOperation<Type>
{
public:
  OperationUpdate(SolutionSpaceBase<Type>& space):
    space(space)
  {};
  ~OperationUpdate(){};
  void Apply(){space.NewtonIterate();};
private:
  SolutionSpaceBase<Type>& space;
};

template <class Type>
class SolutionOrdering
{
public:
  template <class Type2>
  SolutionOrdering<Type>& operator= (const SolutionOrdering<Type2>& orderingToCopy); 

  ~SolutionOrdering();
  UInt size();
  Int Insert(const std::string lineCommand);
  Int Finalize(std::vector<SolutionSpaceBase<Type>*>& solSpaces);
  void Print();
  void Iterate();
  Int ReadOrdering(std::ifstream& fin, std::streampos locBegin, std::streampos locEnd);
  std::vector<std::string> commands;
private:
  std::vector<SolutionOperation<Type>*> opList;
};

template <class Type>
Int ReadSolutionOrdering(SolutionOrdering<Type>& operations, std::string casename, std::string pathname);

//include implementations
#include "solutionOperations.tcc"

#endif
