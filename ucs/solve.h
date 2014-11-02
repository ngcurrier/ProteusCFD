#ifndef SOLVE_H__
#define SOLVE_H__

#include "general.h"
#include "solutionSpaceBase.h"
#include "solutionSpace.h"
#include "solutionOperations.h"
#include "temporalControl.h"
#include <vector>

template <class Type>
void Solve(std::vector<SolutionSpaceBase<Type>*>& solSpaces, SolutionOrdering<Type>& solOperations);

template <class Type>
void ExplicitSolve(SolutionSpace<Type>& solSpace);

enum SolverEnum
  {
    Explicit,
    SGS,
    //GMRES,
    NUMBER_OF_SOLVERS
  };

//include implementations3
#include "solve.tcc"

#endif
