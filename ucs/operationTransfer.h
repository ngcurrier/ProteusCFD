#ifndef OPERATION_TRANSFER__
#define OPERATION_TRANSFER__

// class defines the transfer operation

#include "transfer.h"
#include "fieldTransfer.h"
#include "solutionSpace.h"
#include "solutionOperations.h"


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

//include implementations
#include "operationTransfer.tcc"

#endif
