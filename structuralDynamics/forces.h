#ifndef SFORCES_H__
#define SFORCES_H__

#include <fstream>
#include <iostream>
#include <cmath>

namespace STRUCTDYN{

class SForces
{
public:
  int nnodes;
  int nsurfs;
  int ndof;
  int ndofForces;

  //backward compatibility for 1d case
  bool isOneD;
  double omega;
  double mag;

  int* node;
  int* doftag;
  double* value;

  //memory init flag
  bool init;

  SForces();
  ~SForces();

  void Init(int nnodes, int nsurfs, int dof);
  void Read();
  void Compute(double* forces, double t);
  void Set1DCase(double omega, double mag);

private:
};

}
#endif
