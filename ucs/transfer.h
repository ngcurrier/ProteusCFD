#ifndef TRANSFER_H__
#define TRANSFER_H__

#include "general.h"
#include "solutionSpaceBase.h"

/*
  This class is designed to hold information necessary for transfer from
  one solutionSpaceBase type object to another.  
*/

template <class Type>
class Transfer
{
public:
  Transfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest);
  virtual ~Transfer();
  virtual void DoTransfer() = 0;
private:
  //store base solution space transferring from
  SolutionSpaceBase<Type>* sourceSpace;
  //store base solution space transferring to
  SolutionSpaceBase<Type>* destSpace;
};

template <class Type>
class ScalarTransfer : public Transfer<Type>
{
public:
  ScalarTransfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest, 
		 ScalarField<Type>& sourceField, ScalarField<Type>& destField);
  void DoTransfer();
private:
  //store single value field transferring from
  ScalarField<Type>& sourceScalarField;
  //store single value field transferring to
  ScalarField<Type>& destScalarField;
};


//include implementation
#include "transfer.tcc"

#endif
