#ifndef SPARAM_H__
#define SPARAM_H__

#include "element_lib.h"
#include "forces.h"
#include <iostream>

namespace STRUCTDYN{

class Element;

class SParam
{
public:
  SParam();
  ~SParam();

  int ReadMesh(std::string fileName, SForces* forces);
  void ReadICs();
  void SetICs(double* icx, double* icxd);
  //pass in a displacement vector to move rxyz array
  //returns new vector rxyz of new locations, also can use scale factor
  //to make movements more exaggerated
  void MoveMesh(double* dxyz, double* rxyz, double scale);
  
  //total degrees of freedom
  int dof;

  //number of materials we need to store properties for
  int nMaterials;

  // Flag as to whether the beam rigidity is
  // given, instead of material and cross-sect
  // properties.... If set to TRUE: 
  // Ix = Torsional Rigidity (G*Ix)
  // Iy = Bending Rigidity (E*Iy)
  // Iz = Bending Rigidity (E*Iz)
  // Area = Axial Rigidity (E*A)
  // Rho = Beam Mass (dens*area*length)
  // em = Elastic Modulus = 1.0
  // sm = Shear Modulus = 1.0
  bool* rigidity;

  //ym = Young's Modulus
  double* ym;
  //ar = cross sectional area
  double* ar;
  //gm = shear Modulus of elasticity
  double* sm;
  //rho = density
  double* rho;
  //ix, iy & iz = moment of inertia about x, y and z axis, respectively
  double* ix;
  double* iy;
  double* iz;

  //total number of nodes in mesh
  int nnodes;

  //keep the grid here for now
  double* xyz;

  //keep the elements here for now
  Element* elements;

  //list of element counts
  //0 - beam
  //...
  int nelem[1];

  //check whether to clear memory or not
  bool init;
  //check for memory allocation for ICs
  bool initIC;

  //memory for IC storage
  int ndofICsD;
  int ndofICsV;
  int* icsnode;
  int* icsdoftag;
  double* icsvalue;

  //list nnodes long with calculated dof offsets
  int* nodeOffsetsDOF;
  
private:
};

}

#endif
