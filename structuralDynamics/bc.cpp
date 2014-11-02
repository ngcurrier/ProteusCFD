#include "bc.h"

namespace STRUCTDYN{

BC::BC()
{
  this->init = false;
  this->initstore = false;

  param = NULL;

  return;
}

BC::~BC()
{
  if(this->init){
    delete [] node;
    delete [] doftag;
    delete [] value;
  }
  if(this->initstore){
    delete [] matrixStore;
  }
  return;
}

void BC::SetParam(SParam* param)
{
  this->param = param;
  return;
}

void BC::Read()
{
  int i;
  //currently this only reads in the fixed bcs
  std::string fileName = "infile.bc";
  std::ifstream fin;

  fin.open(fileName.c_str());

  this->init = true;

  fin >> this->ndofbc;

  node = new int[this->ndofbc];
  doftag = new int[this->ndofbc];
  value = new double[this->ndofbc];

  for(i = 0; i < this->ndofbc; i++){
    fin >> node[i];
    fin >> doftag[i];
    fin >> value[i];
  }

  fin.close();

}

void BC::ModifyMatrix(double* keff, int dof)
{
  int i, j, id;

  //============================================
  //
  //
  //WARNING: this is only good for beam elements
  //
  //
  //
  //============================================

  //allocate memory to store the rows we are about to change
  this->matrixStore = new double[dof*ndofbc];
  this->initstore = true;

  //first loop through the list setting fixed values for DOF in list
  //and contributing appropriate pieces to the rhs, also set DOF
  //diagonal to unity 
  for(i = 0; i < ndofbc; i++){
    id = param->nodeOffsetsDOF[node[i]] + doftag[i];
    for(j = 0; j < dof; j++){
      matrixStore[i*dof + j] = keff[j*dof + id];
    }
    //zero the row
    for(j = 0; j < dof; j++){
      keff[id*dof + j] = 0.0;
    }
    //zero the column
    for(j = 0; j < dof; j++){
      keff[j*dof + id] = 0.0;
    }
    //set diagonal to unity
    keff[id*dof + id] = 1.0;
  }


  return;
}

void BC::Apply(double* rhs, int dof)
{
  int i, j, id;

  //============================================
  //
  //
  //WARNING: this is only good for beam elements
  //
  //
  //
  //============================================

  //add contributions to the rhs 
  for(i = 0; i < ndofbc; i++){
    for(j = 0; j < dof; j++){
      rhs[j] -= matrixStore[i*dof + j]*value[i];
    }
    //hard set the rhs to the value we want to enforce
    id = param->nodeOffsetsDOF[node[i]] + doftag[i];
    rhs[id] = value[i];
  }
  
  return;
}

}
