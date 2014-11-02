#ifndef BC_H__
#define BC_H__

#include "structparam.h"
#include <iostream>
#include <fstream>

namespace STRUCTDYN{

class SParam;

class BC
{
public:
  BC();
  ~BC();
  
  //store param object for DOF lookup
  SParam* param;

  //nodes ids for read boundary conditions
  int* node;
  //dof of freedom for node read
  int* doftag;
  //values of specified bc
  double* value;
  
  //store columns of dof which have been modified
  double* matrixStore;

  //memory init flag
  bool init;
  bool initstore;

  //number of nodes with a specified bc in the mesh
  int ndofbc;

  //set parameter to object for internal use
  void SetParam(SParam* param);

  //initialize bc object
  void Read();

  //applies bcs to matrix and rhs given
  void ModifyMatrix(double* keff, int dof);

  //applies bcs to the rhs matrix from stored info
  void Apply(double* rhs, int dof);

private:

};

}


#endif
