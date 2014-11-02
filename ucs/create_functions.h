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

template <class Type>
Int CreateEqnSet(SolutionSpace<Type>* space)
{
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

#endif
