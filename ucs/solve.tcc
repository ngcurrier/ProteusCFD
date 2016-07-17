#include "bc.h"
#include "solutionField.h"
#include "forces.h"
#include "mesh.h"
#include "eqnset.h"
#include "timestep.h"
#include "jacobian.h"
#include "matrix.h"
#include "threaded.h"
#include "parallel.h"
#include "sensors.h"
#include "forces.h"
#include "crs.h"
#include "walldist.h"
#include "timer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

//this function handles implicit solution driving, inner and outer loops, and coupling
template <class Type>
void Solve(std::vector<SolutionSpaceBase<Type>*>& solSpaces, SolutionOrdering<Type>& solutionOperations)
{
  TemporalControl<Real>& temporalControl = solSpaces[0]->temporalControl;

  //this loop handles things that are done only once per solver run before timestepping
  for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
      it != solSpaces.end(); ++it){
    SolutionSpaceBase<Type> & space = **it;
    space.PreIterate();
  }

  std::cout << "Beginning core flow solver:" << std::endl;
  std::cout << "===========================\n" << std::endl;

  //Outer timestepping loop
  Int iter = 1;
  while(iter <= temporalControl.nSteps){

    //Do things that should happen pre timestep
    for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
	it != solSpaces.end(); ++it){
      SolutionSpaceBase<Type> & space = **it;
      space.PreTimeAdvance();
    }

    //Do the Newton inner iteration
    for(Int newtonit = 1; newtonit <= temporalControl.newtonIter; newtonit++){
      //this iterates over things defined in the param file in order
      //such as advancing solution, transfers, etc.
      solutionOperations.Iterate();
    }
    
    //Do things that should happen post timestep
    for(typename std::vector<SolutionSpaceBase<Type>*>::iterator it = solSpaces.begin(); 
	it != solSpaces.end(); ++it){
      SolutionSpaceBase<Type> & space = **it;
      space.PostTimeAdvance();
    }
    
    iter++;
  }

  return;
}

//this function handles explicit solution methodology only and contains special
//logic for dealing with solution spaces that do not naturally use conservative
//variables for storage/update - i.e. those that are preconditioned
template <class Type>
void ExplicitSolve(SolutionSpace<Type>& solSpace)
{
  Int i, j;
  EqnSet<Type>* eqnset = solSpace.eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = eqnset->nauxvars + neqn;
  Type cn = 1.0;
  Mesh<Type>* m = solSpace.m;
  CRS<Type>& crs = *solSpace.crs;
  Int nnode = m->GetNumNodes();

  Type* dt = solSpace.GetField("timestep", FIELDS::STATE_NONE);

#ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(neqn, m, eqnset, dt, cn, nvars) \
  private(i, j) 
#endif
    {
      if(eqnset->varsConservative){
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
	for(i = 0; i < nnode; i++){
	  for(j = 0; j < neqn; j++){
	    crs.x[i*neqn + j] = crs.b[i*neqn + j]*dt[i]/m->vol[i];
	  }
	  eqnset->ApplyDQ(&crs.x[i*neqn], &solSpace.q[i*nvars], &m->xyz[i]);
	}
      }
      else{
#if 0
	Type* betaa = solSpace.GetField("beta", FIELDS::STATE_NONE);
	Type* A = (Type*)alloca(sizeof(Type)*neqn*neqn);
	Int p[neqn];
	//This works but is not currently functional for any boundary condition where the 
	//jacobian and residual must be modified to enforce the BC
	for(i = 0; i < nnode; i++){
	  Type beta = betaa[i];
	  MemBlank(A, neqn*neqn);
	  eqnset->ContributeTemporalTerms(&solSpace.q[i*nvars], m->vol[i], cn, dt[i], 0.0, A, beta);
	  LU(A, p, neqn);
	  LuSolve(A, &crs.b[i*neqn], p, &crs.x[i*neqn], neqn);
	  //solution is return in res.. copy to correct place
	  memcpy(&crs.x[i*neqn], &crs.b[i*neqn], neqn*sizeof(Type));
	}
#else
	Type* qc = new Type[neqn];
	for(i = 0; i < nnode; i++){
	  memcpy(qc, &solSpace.q[i*nvars], sizeof(Type)*neqn);
	  eqnset->NativeToConservative(&solSpace.q[i*nvars]);
	  for(j = 0; j < neqn; j++){
	    crs.x[i*neqn + j] = crs.b[i*neqn + j]*dt[i]/m->vol[i];
	  }
	  for(j = 0; j < neqn; j++){
	    solSpace.q[i*nvars + j] += crs.x[i*neqn + j];
	  }
	  eqnset->ConservativeToNative(&solSpace.q[i*nvars]);
	  //get the change in the native variable type
	  for(j = 0; j < neqn; j++){
	    crs.x[i*neqn + j] = solSpace.q[i*nvars + j] - qc[j];
	  }
	  memcpy(&solSpace.q[i*nvars], qc, sizeof(Type)*neqn);
	}
	delete [] qc;
#endif
      }
    }
    return;
}
