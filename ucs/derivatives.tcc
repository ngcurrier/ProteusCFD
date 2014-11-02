#include "forces.h"
#include "mesh.h"  
#include "eqnset.h"
#include "sensors.h"
#include "portFileio.h"
#include "exceptions.h"

template <class Type>
Type Compute_Obj_Function(SolutionSpace<Type>& space)
{
  Int i;
  Mesh<Type>* m = space.m;
  EqnSet<Type>* eqnset = space.eqnset;
  Forces<Type>* forces = space.forces;
  Param<Type>* param = space.param;
  space.forces->Compute();
  space.sensors->Search();
  space.sensors->Sample();

  Type obj;

  if(param->objFuncId == 0){
    obj = 1.0/forces->bodies[1].cl;
  }
  else if(param->objFuncId == 1){
    // lift/drag
    //
    Type frac = forces->bodies[1].cd/forces->bodies[1].cl;
    obj = 0.5*frac*frac;
  }
  else if(param->objFuncId == 2){
    // lift/drag + moment
    //
    Type frac = CAbs(0.7*forces->bodies[1].cd/forces->bodies[1].cl) + 0.3*CAbs(forces->bodies[1].cm);
    obj = 0.5*frac*frac;
  }
  else if(param->objFuncId == 3){
    // this is a least squares sensor target objective function
    //
    Type LSQ;
    Type diff;
    Sensors<Type>* sens = space.sensors;
    std::vector<Type>& sensTarget = param->sensTarget;
    Int nsensors = sensTarget.size();
    if(nsensors != sens->count){
      Abort << "WARNING: Number of sensors and number of targets are mismatched";
    }
    Int neqn = eqnset->neqn;			\
    Int eqnseek = param->sensEqn;

    LSQ = 0.0;
    for(i = 0; i < nsensors; i++){
      //use a relative difference otherwise small changes at low concentration sensors
      //get steamrolled
      Type value = sensTarget.at(i);
      //diff = (sens->values[i*neqn + eqnseek] - value)/value;
      diff = sens->values[i*neqn + eqnseek] - value;
      diff = diff*diff;
      LSQ += diff;
    }
    LSQ = sqrt(LSQ);
    obj = LSQ;
  }
  else{
    Abort << "WARNING: Objective function not defined";
  }

  return (obj);
}
