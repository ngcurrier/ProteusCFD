
template <class Type>
Transfer<Type>::Transfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest) :
  sourceSpace(source), destSpace(dest)
{}

template <class Type>
Transfer<Type>::~Transfer()
{}

template <class Type>
ScalarTransfer<Type>::ScalarTransfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest, 
				     ScalarField<Type>& sourceField, ScalarField<Type>& destField) :
  Transfer<Type>(source, dest), sourceScalarField(sourceField), destScalarField(destField) 
{}
  
template <class Type>
void ScalarTransfer<Type>::DoTransfer()
{
  destScalarField.SetField(sourceScalarField.GetField());
}

