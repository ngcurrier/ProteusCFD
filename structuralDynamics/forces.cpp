#include "forces.h"

namespace STRUCTDYN{

SForces::SForces()
{
  this->init = false;
  return;
}

SForces::~SForces()
{
  if(this->init){
    delete [] node;
    delete [] doftag;
    delete [] value;
  }
  return;
}

void SForces::Init(int nnodes, int nsurfs, int ndof)
{
  this->nnodes = nnodes;
  this->nsurfs = nsurfs;
  this->ndof = ndof;

  isOneD = false;

  return;
}

void SForces::Read()
{
  int i;
  //currently this only reads in the fixed bcs
  std::string fileName = "infile.forces";
  std::ifstream fin;

  fin.open(fileName.c_str());

  this->init = true;

  fin >> this->ndofForces;

  node = new int[this->ndofForces];
  doftag = new int[this->ndofForces];
  value = new double[this->ndofForces];

  for(i = 0; i < this->ndofForces; i++){
    fin >> node[i];
    fin >> doftag[i];
    fin >> value[i];
  }

  fin.close();


  return;
}

void SForces::Set1DCase(double omega, double mag)
{
  this->omega = omega;
  this->mag = mag;
  isOneD = true;

  return;
}

void SForces::Compute(double* forces, double t)
{
  int i, id;
  int dofbeam = 6;

  //assume zero forces unless otherwise noted in forces file
  for(i = 0; i < ndof; i++){
    forces[i] = 0.0;
  }
  
  //
  //WARNING: this only works for 3D beam elements right now
  //
  for(i = 0; i < ndofForces; i++){
    id = node[i]*dofbeam + doftag[i];
    forces[id] = value[i];
  }

  if(isOneD){
    double phip = 0.0;
    for(i = 0; i < ndof; i++){
      forces[i] = mag*sin(omega*t - phip);
    }
  }

  return;
}

}
