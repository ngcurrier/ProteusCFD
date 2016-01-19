#include "structparam.h"

namespace STRUCTDYN{

SParam::SParam()  
{
  this->init = false;
  this->initIC = false;
  //dumb allocator
  return;
}

SParam::~SParam()
{
  //free the memory
  if(this->init){
    delete [] rigidity;
    delete [] ym;
    delete [] ar;
    delete [] sm;
    delete [] rho;
    delete [] ix;
    delete [] iy;
    delete [] iz;
    delete [] xyz;
    delete [] elements;
    delete [] nodeOffsetsDOF;
  }
  if(this->initIC){
    delete[] icsnode;
    delete[] icsdoftag;
    delete[] icsvalue;
  }

  return;
}

int SParam::ReadMesh(std::string fileName, SForces* forces)
{
  int i, j;
  int count;
  int nCfd;
  int nid;
  int elemId;
  int matId;

  std::ifstream fin;
  fin.open(fileName.c_str());
  if(!fin.is_open()){
    std::cerr << "Could not open mesh file " << fileName << std::endl;
    return 1;
  }
  
  //read number of nodes
  fin >> this->nnodes;

  this->nodeOffsetsDOF = new int[this->nnodes];

  //read in number of materials we have to track
  fin >> nMaterials;

  //read number of cfd nodes on surface
  fin >> nCfd;

  //read number of beam elements
  fin >> this->nelem[0];
  
  //create elements
  this->elements = new Beam[nelem[0]]();

  //allocate coordinates memory
  xyz = new double[nnodes*3];

  //read xyz coordinates of the nodes
  for(i = 0; i < nnodes; i++){
    //read in node id
    fin >> nid;
    for(j = 0; j < 3; j++){
      fin >> xyz[nid*3 + j];
    }
  }
  
  //read element connectivity
  int nodes[2];
  for(i = 0; i < nelem[0]; i++){
    int matId;
    double vecxy[3];
    fin >> elemId;
    int nelemnodes = this->elements[elemId].GetNnodes();
    for(j = 0; j < nelemnodes; j++){
      fin >> nodes[j];
    }
    fin >> matId;
    for(j = 0; j < 3; j++){
      fin >> vecxy[j];
    }
    elements[elemId].Init(matId, vecxy, nodes);
  }

  //initialize node offset list
  for(i = 0; i < nnodes; i++){
    nodeOffsetsDOF[i] = -1;
  }

  //setup the offset list, simply count over all the elements
  this->dof = 0;
  count = 0;
  for(i = 0; i < nelem[0]; i++){
    int ndofpn = elements[i].GetNodeDOF();
    for(j = 0; j < elements[i].GetNnodes(); j++){
      if(nodeOffsetsDOF[elements[i].nodes[j]] >= 0){
	//do nothing we've already added in this node
      }
      else{
	//add the dof for the node
	nodeOffsetsDOF[elements[i].nodes[j]] = count;
	count += ndofpn;
	this->dof += ndofpn;
      }
    }
  }

  //allocate the space for storing material properties
  rigidity = new bool[nMaterials];
  ym = new double[nMaterials];
  ar = new double[nMaterials];
  sm = new double[nMaterials];;
  rho = new double[nMaterials];
  ix = new double[nMaterials];
  iy = new double[nMaterials];
  iz = new double[nMaterials];

  //read in material properties
  for(i = 0; i < nMaterials; i++){
    fin >> matId;
    fin >> rigidity[matId];
    fin >> ym[matId];
    fin >> sm[matId];
    fin >> rho[matId];
    fin >> ar[matId];
    fin >> ix[matId];
    fin >> iy[matId];
    fin >> iz[matId];
  }
  
  this->init = true;

  //initialize forces object
  int nsurf = 1;
  forces->Init(nnodes, nsurf, dof);

  fin.close();

  return dof;
}

void SParam::ReadICs()
{
  int i;
  //currently this only reads in the fixed bcs
  std::string fileName = "infile.ics";
  std::ifstream fin;

  fin.open(fileName.c_str());

  this->initIC = true;

  //number of ICs on displacement
  fin >> this->ndofICsD;
  //number of ICs on velocity
  fin >> this->ndofICsV;

  icsnode = new int[this->ndofICsD+this->ndofICsV];
  icsdoftag = new int[this->ndofICsD+this->ndofICsV];
  icsvalue = new double[this->ndofICsD+this->ndofICsV];

  for(i = 0; i < this->ndofICsD+this->ndofICsV; i++){
    fin >> icsnode[i];
    fin >> icsdoftag[i];
    fin >> icsvalue[i];
  }

  fin.close();


  return;
}

void SParam::SetICs(double* icx, double* icxd)
{
  int i, id;
  int dofbeam = 6;

  ReadICs();
  
  //assume ICs are zero unless otherwise given
  for(i = 0; i < dof; i++){
    icx[i] = icxd[i] = 0.0;
  }
  
  //
  //WARNING: this is only valid for 3D beam elements at this time
  //
  for(i = 0; i < ndofICsD; i++){
    id = icsnode[i]*dofbeam + icsdoftag[i];
    icx[id] = icsvalue[i];
  }
  for(i = ndofICsV; i < ndofICsD+ndofICsV; i++){
    id = icsnode[i]*dofbeam + icsdoftag[i];
    icxd[id] = icsvalue[i];
  }
	
  return;
}

void SParam::MoveMesh(double* dxyz, double* rxyz, double scale)
{
  int i, j;
  for(i = 0; i < nnodes; i++){
    for(j = 0; j < 3; j++){
      rxyz[i*3 + j] = scale*dxyz[i*3 + j] + xyz[i*3 + j];
    }
  }
  return;
}

}
