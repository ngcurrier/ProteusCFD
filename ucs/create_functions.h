#ifndef CREATE_FUNCTIONS_H__
#define CREATE_FUNCTIONS_H__

#include "general.h"
#include "eqnset.h"
#include "compressible.h"
#include "incompressible.h"
#include "compressibleFR.h"
#include "heatTransfer.h"
#include "mesh.h"
#include "param.h"
#include "exceptions.h"
#include "solutionSpace.h"
#include "fluid_structure.h"
#include <vector>

template <class Type>
Int CreateEqnSet(SolutionSpace<Type>* space)
{
  std::cout << "CREATING COPY OF EQUATION SET CreateEqnSet(SolutionSpace<>* space)" << std::endl;
  if(space == NULL){
    std::cerr << "Solution Space point is NULL in CreateEqnSet()... FAILING!" << std::endl;
    return (-1);
  }

  //setup correct eqnset object..call also allocates solution memory
  if(space->param->eqnset_id == CompressibleEuler || space->param->eqnset_id == CompressibleNS){
    space->eqnset = new CompressibleEqnSet<Type>(space, space->param);
  }
  else if(space->param->eqnset_id == IncompressibleEuler || space->param->eqnset_id == IncompressibleNS){
    space->eqnset = new IncompressibleEqnSet<Type>(space, space->param);
  }
  else if(space->param->eqnset_id == CompressibleEulerFR || space->param->eqnset_id == CompressibleNSFR){
    space->eqnset = new CompressibleFREqnSet<Type>(space, space->param);
  }
  else if(space->param->eqnset_id == HeatTransfer){
    space->eqnset = new HeatTransferEqnSet<Type>(space, space->param);
  }
  else{
    std::cerr << "Solver type " << space->param->eqnset_id << " not found.. FATAL\n";
    Abort << "Dead";
    return (-1);
  }
  return (0);
};


//this is simply a utility function which gives us access to a semi-functional
//eqnset object for things like complex jacobians
template <class Type>
Int CreateEqnSet(EqnSet<Type>** eqnset, Param<Type>* param){
  std::cout << "CREATING DUMB COPY OF EQUATION SET CreateEqnSet(EqnSet<>**, Param<>* param)" << std::endl;

  //setup correct eqnset object..call also allocates solution memory
  if(param->eqnset_id == CompressibleEuler || param->eqnset_id == CompressibleNS){
    *eqnset = new CompressibleEqnSet<Type>(NULL, param);
  }
  else if(param->eqnset_id == IncompressibleEuler || param->eqnset_id == IncompressibleNS){
    *eqnset = new IncompressibleEqnSet<Type>(NULL, param);
  }
  else if(param->eqnset_id == CompressibleEulerFR || param->eqnset_id == CompressibleNSFR){
    *eqnset = new CompressibleFREqnSet<Type>(NULL, param);
  }
  else if(param->eqnset_id == HeatTransfer){
    *eqnset = new HeatTransferEqnSet<Type>(NULL, param);
  }
  else{
    std::cerr << "Solver type " << param->eqnset_id << " not found.. FATAL!\n";
    Abort << "Dead";
    return (-1);
  }
  return (0);
}

// This is the factory for solution spaces from param files read
template <class Type>
std::vector<SolutionSpaceBase<Type>*> CreateSolutionSpaces(typename std::vector<Param<Type>* > paramList, PObj<Type>& pobj,
							   TemporalControl<Type>& temporalControl){
  std::vector<SolutionSpaceBase<Type>*> solSpaces;
  for(typename std::vector<Param<Type>* >::iterator it = paramList.begin();
      it != paramList.end(); ++it){
    Param<Type>* param = *it;
    SolutionSpaceBase<Type>* solSpace;
    if(param->spacename == "structure"){
      solSpace = new STRUCTDYN::Structure(param, param->spacename, temporalControl);
    }
    else if(param->spacename == "typicalSection"){
      solSpace = new STRUCTDYN::TypicalSection(param, param->spacename, temporalControl);
    }
    else{
      solSpace = new SolutionSpace<Type>(param, &pobj, param->spacename, temporalControl);
    }
    solSpaces.push_back(solSpace);
  }
  for(typename std::vector<SolutionSpaceBase<Type>* >::iterator it = solSpaces.begin(); 
      it != solSpaces.end(); ++it){
    SolutionSpaceBase<Type> & space = **it;
    space.WriteAvailableFields();
  }
  return solSpaces;
}

#endif
