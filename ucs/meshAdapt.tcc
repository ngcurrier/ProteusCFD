template <class Type>
void Mesh<Type>::RefineNodeBased(Bool* refineNodeList)
{
  return;
}

//this routine must resize all arrays in MemInitMesh()
template <class Type>
void Mesh<Type>::RebuildRefinedMesh(Type* xyzNew, Int xyzNewCount, 
				    Element<Type>** newElements, Int newElemCount, 
				    Int* deadElemList, Int deadElemCount)
{
  return;
}
