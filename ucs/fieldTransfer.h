#ifndef FIELD_TRANSFER_H__
#define FIELD_TRANSFER_H__

#include "general.h"
#include "transfer.h"
#include "solutionSpaceBase.h"
#include "solutionField.h"


template <class Type>
class FieldTransfer : public Transfer<Type>
{
public:
  FieldTransfer(SolutionSpaceBase<Type>* source, SolutionSpaceBase<Type>* dest, 
		SolutionField<Type>& sourceField, SolutionField<Type>& destField);
  ~FieldTransfer(){};
  void DoTransfer();
  void SetXyz(Type* sourceXyz, Int nptsSource, Type* destXyz, Int nptsDest);
  //  this is optional to set but will increase accuracy, locality and speed
  //  ipsp - (crs indx into point surrounding point map); psp - (point to point map)
  void SetPointToPointMap(Int* ipspSource, Int* pspSource);
private:
  //store array field on a mesh transferring from
  SolutionField<Type>& sourceArrayField;
  //store coords. transferring from, relates 1 to 1 with array field
  Type* sourceXyz;
  Int nptsSource;

  //store array field on a mesh transferring to
  SolutionField<Type>& destArrayField;
  //store coords. transferring from, relates 1 to 1 with array field
  Type* destXyz;
  Int nptsDest;

  //Point surrounding point CRS maps
  Int* ipspSource;
  Int* pspSource;

  Bool geometrySet;
};

//include implementation
#include "fieldTransfer.tcc"

#endif
