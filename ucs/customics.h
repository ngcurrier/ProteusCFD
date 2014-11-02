#ifndef CUSTOM_ICS_H__
#define CUSTOM_ICS_H__

#include <iostream>
#include "general.h"
#include "mesh.h"
#include "param.h"
#include "solutionSpace.h"

template <class Type>
void InitIC(Int icid, SolutionSpace<Type>* space);

//include implementations
#include "customics.tcc"

#endif
