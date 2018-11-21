
template <class Type>
OperationTransfer<Type>::OperationTransfer(SolutionSpaceBase<Type>& source, SolutionSpaceBase<Type>& dest, 
					   ScalarField<Type>& fieldSource, ScalarField<Type>& fieldDest)
{
  //create a scalar transfer
  transfer = new ScalarTransfer<Type>(&source, &dest, fieldSource, fieldDest);
}

template <class Type>
OperationTransfer<Type>::OperationTransfer(SolutionSpaceBase<Type>& source, SolutionSpaceBase<Type>& dest, 
					   SolutionField<Type>& fieldSource, SolutionField<Type>& fieldDest)
{
  //create a field transfer
  transfer = new FieldTransfer<Type>(&source, &dest, fieldSource, fieldDest);
}

template <class Type>
OperationTransfer<Type>::~OperationTransfer()
{
  delete transfer;
}

template <class Type>
void OperationTransfer<Type>::Apply()
{
  transfer->DoTransfer();
}

