#include "element_lib.h"

namespace STRUCTDYN{

void Element::Transform(double* elm, double* eldiagm, double* els, bool dynamic, SParam* param)
{
  int i, j;
  int ndofpe = this->GetElemDOF();
  double* tran = new double[ndofpe*ndofpe];
  double* trant = new double[ndofpe*ndofpe];
  double* temp = new double[ndofpe*ndofpe];

  //zero transform matrices
  for(i = 0; i < ndofpe*ndofpe; i++){
    tran[i] = trant[i] = 0.0;
  }

  //compute global coordinate vectors
  int n1 = this->nodes[0];
  int n2 = this->nodes[1];

  double xhat[3];
  double yhat[3];
  double zhat[3];
  double aux[3];
  
  double* xyz = param->xyz;
  xhat[0] = xyz[n2*3 + 0] - xyz[n1*3 + 0];
  xhat[1] = xyz[n2*3 + 1] - xyz[n1*3 + 1];
  xhat[2] = xyz[n2*3 + 2] - xyz[n1*3 + 2];

  aux[0] = vecxy[0];
  aux[1] = vecxy[1];
  aux[2] = vecxy[2];

  double d;
  d = Normalize(xhat, xhat);
  d = Normalize(aux, aux);

  //compute cross product of xhat and x-y plane vector to get zhat
  CrossProduct(xhat, aux, zhat);
  d = Normalize(zhat, zhat);
  
  //comptue cross product of zhat and xhat to get yhat
  CrossProduct(zhat, xhat, yhat);
  d = Normalize(yhat, yhat);

  //form the transformation matrix
  for(j = 0; j < 3; j++){
    tran[0*ndofpe + j] = xhat[j];
    tran[1*ndofpe + j] = yhat[j];
    tran[2*ndofpe + j] = zhat[j];
  }
  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++){
      tran[(i+3)*ndofpe + (j+3)] = tran[i*ndofpe + j];
      tran[(i+6)*ndofpe + (j+6)] = tran[i*ndofpe + j];
      tran[(i+9)*ndofpe + (j+9)] = tran[i*ndofpe + j];
    }
  }
  //compute transpose
  for(i = 0; i < ndofpe; i++){
    for(j = 0; j < ndofpe; j++){
      trant[j*ndofpe + i] = tran[i*ndofpe + j];
    }
  }

  //perform triple product to put element mass and stiffness matrices
  //in global coordinates
  MatMatMult(trant, els, temp, ndofpe);
  MatMatMult(temp, tran, els, ndofpe);

  if(dynamic){
    MatMatMult(trant, elm, temp, ndofpe);
    MatMatMult(temp, tran, elm, ndofpe);
    if(eldiagm != NULL){
      MatMatMult(trant, eldiagm, temp, ndofpe);
      MatMatMult(temp, tran, eldiagm, ndofpe);
    }
  }


  delete [] tran;
  delete [] trant;
  delete [] temp;

  return;
}

Element::Element()
{
  //don't call this
  return;
}

Element::~Element()
{
  delete [] nodes;
  return;
}

void Element::Init(int materialId, double* vecxy, int* nodes)
{
  int i;

  //copy over constructor arguments
  this->materialId = materialId;


  for(i = 0; i < 3; i++){
    this->vecxy[i] = vecxy[i];
  }

  this->nodes = new int[GetNnodes()];
  for(i = 0; i < GetNnodes(); i++){
    this->nodes[i] = nodes[i];
  }
  
  return;
}

Beam::Beam()
{
  return;
}


void Beam::Compute(double* elm, double* els, bool dynamic, bool gravity, SParam* param)
{
  int i, j;
  
  int ndofpe = GetElemDOF();

  double ix = param->ix[this->materialId];
  double iy = param->iy[this->materialId];
  double iz = param->iz[this->materialId];
  
  double ar = param->ar[this->materialId];
  double ym = param->ym[this->materialId];
  double sm = param->sm[this->materialId];
  double den = param->rho[this->materialId];
  
  int rigidity = param->rigidity[this->materialId];

  if(rigidity){
    sm = 1.0;
    ym = 1.0;
  }

  double* xyz = param->xyz;

  int n1 = this->nodes[0];
  int n2 = this->nodes[1];
  double dx = xyz[n1*3 + 0] - xyz[n2*3 + 0];
  double dy = xyz[n1*3 + 1] - xyz[n2*3 + 1];
  double dz = xyz[n1*3 + 2] - xyz[n2*3 + 2];
  
  double xl = sqrt(dx*dx + dy*dy + dz*dz);

  //zero matrices
  for(i = 0; i < ndofpe; i++){
    for(j = 0; j < ndofpe; j++){
      elm[i*ndofpe + j] = 0.0;
      els[i*ndofpe + j] = 0.0;
    }
  }

  //compute element stiffness matrix
  int n = GetElemDOF();
  els[0*n + 0] = ym*ar/xl;
  els[0*n + 6] = -els[0*n + 0];
  els[1*n + 1] = 12.0*ym*iz/(xl*xl*xl);
  els[1*n + 5] = 6.0*ym*iz/(xl*xl);
  els[1*n + 7] = -els[1*n + 1];
  els[1*n + 11] = els[1*n + 5];
  els[2*n + 2] = 12.0*ym*iy/(xl*xl*xl);
  els[2*n + 4] = -6.0*ym*iy/(xl*xl);
  els[2*n + 8] = -els[2*n + 2];
  els[2*n + 10] = els[2*n + 4];
  els[3*n + 3] = (sm*ix)/xl;
  els[3*n + 9] = -els[3*n + 3];
  els[4*n + 4] = 4.0*ym*iy/xl;
  els[4*n + 8] = 6.0*ym*iy/(xl*xl);
  els[4*n + 10] = 2.0*ym*iy/xl;
  els[5*n + 5] = 4.0*ym*iz/xl;
  els[5*n + 7] = -6.0*ym*iz/(xl*xl);
  els[5*n + 11] = 2.0*ym*iz/xl;
  els[6*n + 6] = ym*ar/xl;
  els[7*n + 7] = 12.0*ym*iz/(xl*xl*xl);
  els[7*n + 11] = -6.0*ym*iz/(xl*xl);
  els[8*n + 8] = 12.0*ym*iy/(xl*xl*xl);
  els[8*n + 10] = 6.0*ym*iy/(xl*xl);
  els[9*n + 9] = (sm*ix)/xl;
  els[10*n + 10] = 4.0*ym*iy/xl;
  els[11*n + 11] = 4.0*ym*iz/xl;

  //symmetric
  for(i = 0; i < ndofpe; i++){
    for(j = i; j < ndofpe; j++){
      els[j*ndofpe + i] = els[i*ndofpe + j];
    }
  }
  
  double rotf1, rotf2, rotf3, rotf4, rotf5, ral;

  //compute rotatory inertia terms
  if(dynamic || gravity){
    if(rigidity){
      //...only rigidity is given, not intertias nor area. hence
      //   can not compute rotatory inertia. 
      rotf1 = 0.0;
      rotf2 = 0.0;
      rotf3 = 0.0;
      rotf4 = 0.0;
      rotf5 = 0.0;
      ral = den;
    }
    else{
      rotf1 = 6.0/(5.0*ar*xl*xl);
      rotf2 = 1.0/(10.0*ar*xl);
      rotf3 = 2.0/(15.0*ar);
      rotf4 = 1.0/(30.0*ar);
      rotf5 = 1.0/(3.0*ar);
      ral = den*ar*xl;
    }
  }

  //compute element mass matrix
  //consistent mass matrix = mass inertia + rotatory inertia
  //                         (lower freq)   (higher freq)
  elm[0*n + 0] = 1.0/3.0;
  elm[0*n + 6] = 1.0/6.0;
  elm[1*n + 1] = 13.0/35.0 + iz*rotf1;
  elm[1*n + 5] = 11.0*xl/210.0 + iz*rotf2;
  elm[1*n + 7] = 9.0/70.0 - iz*rotf1;
  elm[1*n + 11] = -13.0*xl/420.0 + iz*rotf2;
  elm[2*n + 2] = 13.0/35.0 + iy*rotf1;
  elm[2*n + 4] = -11.0*xl/210.0 - iy*rotf2;
  elm[2*n + 8] = 9.0/70.0 - iy*rotf1;
  elm[2*n + 10] = 13.0*xl/420.0 - iy*rotf2;
  elm[3*n + 3] = ix*rotf5;
  elm[3*n + 9] = 0.50*ix*rotf5;
  elm[4*n + 4] = xl*xl/105.0 + iy*rotf3;
  elm[4*n + 8] = -13.0*xl/420.0 + iy*rotf2;
  elm[4*n + 10] = -xl*xl/140.0 - iy*rotf4;
  elm[5*n + 5] = xl*xl/105.0 + iz*rotf3;
  elm[5*n + 7] = 13.0*xl/420.0 - iz*rotf2;
  elm[5*n + 11] = -xl*xl/140.0 - iz*rotf4;
  elm[6*n + 6] = 1.0/3.0;
  elm[7*n + 7] = 13.0/35.0 + iz*rotf1;
  elm[7*n + 11] = -11.0*xl/210.0 - iz*rotf2;
  elm[8*n + 8] = 13.0/35.0 + iy*rotf1;
  elm[8*n + 10] = 11.0*xl/210.0 + iy*rotf2;
  elm[9*n + 9] = ix*rotf5;
  elm[10*n + 10] = xl*xl/105.0 + iy*rotf3;
  elm[11*n + 11] = xl*xl/105.0 + iz*rotf3;

  //symmetric
  for(i = 0; i < ndofpe; i++){
    for(j = i; j < ndofpe; j++){
      elm[j*ndofpe + i] = elm[i*ndofpe + j];
    }
  }
  
  //multiply by density term
  for(i = 0; i < ndofpe; i++){
    for(j = 0; j < ndofpe; j++){
      elm[i*ndofpe + j] = elm[i*ndofpe + j]*ral;
    }
  }

  return;
}

void Beam::Assemble(double* elm, double* eldiagm, double* els, double* gm, 
		    double* diagm, double* gs, double* gc, bool damping, bool gravity, 
		    double alpha, double beta, SParam* param)
{
  /////////////////////////////////////////////////
  //WARNING!!!!
  //at this point we assume that we have a mesh
  //of all the same element type... this changes if we have
  //mixed elements
  ///////////////////////////////////////////////////

  int i, j, ii;
  int irow, icol;
  int ndofpe = GetElemDOF();
  int ndofpn = GetNodeDOF();
  int id[ndofpe];
  int md[GetNnodes()];
  int tdof = param->dof;


  //find the lower numbered node
  if(nodes[0] < nodes[1]){
    md[0] = nodes[0];
    md[1] = nodes[1];
  }
  else{
    md[0] = nodes[1];
    md[1] = nodes[0];
  }

  //load up array with global id offsets
  ii = 0; 
  for(i = 0; i < GetNnodes(); i++){
    for(j = 0; j < ndofpn; j++){
      id[ii] = param->nodeOffsetsDOF[md[i]] + j;
      ii++;
    }
  }

  for(i = 0; i < ndofpe; i++){
    //fill out the full matrix now... not just the upper part
    for(j = 0; j < ndofpe; j++){
      irow = id[i]*tdof;
      icol = id[j];
      gs[irow + icol] += els[i*ndofpe + j];
      gm[irow + icol] += elm[i*ndofpe + j];
      gc[irow + icol] += alpha*elm[i*ndofpe + j] + beta*els[i*ndofpe + j];
    }
    if(eldiagm != NULL && diagm != NULL){
      for(j = 0; j < ndofpe; j++){
	irow = id[i]*tdof;
	icol = id[j];
	diagm[irow + icol] += eldiagm[i*ndofpe + j];
      }
    }
  }

  return;
}

void Beam::LumpHRZ(double* elm, SParam* param)
{
  int i;

  double ix = param->ix[this->materialId];
  double iy = param->iy[this->materialId];
  double iz = param->iz[this->materialId];
  
  double ar = param->ar[this->materialId];
  double ym = param->ym[this->materialId];
  double sm = param->sm[this->materialId];
  double den = param->rho[this->materialId];
  
  int rigidity = param->rigidity[this->materialId];

  if(rigidity){
    sm = 1.0;
    ym = 1.0;
  }

  double* xyz = param->xyz;

  int n1 = this->nodes[0];
  int n2 = this->nodes[1];
  double dx = xyz[n1*3 + 0] - xyz[n2*3 + 0];
  double dy = xyz[n1*3 + 1] - xyz[n2*3 + 1];
  double dz = xyz[n1*3 + 2] - xyz[n2*3 + 2];
  
  double xl = sqrt(dx*dx + dy*dy + dz*dz);
  
  //number of DOF per element
  int ndofpe = GetElemDOF();
  int n = GetElemDOF();

  double ral = den*ar*xl;
  double rotf1 = 6.0/(5.0*ar*xl*xl);
  double rotf2 = 1.0/(10.0*ar*xl);
  double rotf3 = 2.0/(15.0*ar);
  double rotf4 = 1.0/(30.0*ar);
  double sfac;

  //zero matrix
  for(i = 0; i < ndofpe*ndofpe; i++){
    elm[i] = 0.0;
  }

  //Consistent Mass Matrix = mass inertia + rotatory inertia
  //                         (lower freq)    (higher freq)
  elm[0*n + 0] = 1.0/3.0;
  elm[1*n + 1] = 13.0/35.0 + iz*rotf1;
  elm[2*n + 2] = 13.0/35.0 + iy*rotf1;
  elm[3*n + 3] = ix/(3.0*ar);
  elm[4*n + 4] = xl*xl/105.0 + iy*rotf3;
  elm[5*n + 5] = xl*xl/105.0 + iz*rotf3;
  elm[6*n + 6] = 1.0/3.0;
  elm[7*n + 7] = 13.0/35.0 + iz*rotf1;
  elm[8*n + 8] = 13.0/35.0 + iy*rotf1;
  elm[9*n + 9] = ix/(3.0*ar);
  elm[10*n + 10] = xl*xl/105.0 + iy*rotf3;
  elm[11*n + 11] = xl*xl/105.0 + iz*rotf3;

  for(i = 0; i < ndofpe; i++){
    elm[i*ndofpe + i] = ral*elm[i*ndofpe + i];
  }

  // compute scaling from translational DOFs
  sfac = elm[0*n + 0] + elm[1*n + 1] + elm[2*n + 2] + elm[6*n + 6] + elm[7*n + 7] + elm[8*n + 8];
  sfac = 1.0/sfac;

  for(i = 0; i < ndofpe; i++){
    elm[i*ndofpe + i] = sfac*elm[i*ndofpe + i];
  }

  return;
}

void Beam::InterpolatePoint(double* xyzPt, SParam* param, double* sol, int* ndof,
			    double* values)
{
  int i, j, node;
  double* xyz = param->xyz;
  int nnodes = this->GetNnodes();
  int ndofpn = this->GetNodeDOF();
  *ndof = ndofpn;
  double d;
  //weights
  double* w = new double[nnodes];
  double wsum = 0.0;
  int nodeDOF;
  
  //compute weights
  for(i = 0; i < nnodes; i++){
    d = 0.0;
    for(j = 0; j < 3; j++){
      node = this->nodes[i];
      d = (xyz[node*3 + j] - xyzPt[j])*(xyz[node*3 + j] - xyzPt[j]);
    }
    d = sqrt(d);
    w[i] = 1.0/d;
    wsum += w[i];
  }

  //Normalize weights
  for(i = 0; i < nnodes; i++){
    w[i] = w[i]/wsum;
  }

  //use inverse distance weighting to compute values at new point
  for(i = 0; i < ndofpn; i++){
    values[i] = 0.0;
  }
  for(j = 0; j < nnodes; j++){
    node = this->nodes[j];
    nodeDOF = param->nodeOffsetsDOF[node];
    for(i = 0; i < ndofpn; i++){
      values[i] += sol[nodeDOF + i]*w[j];
    }
  }

  delete [] w;

  return;
}


void Beam::InterpolateValuesToNodes(double* xyzPt, SParam* param, double* values, 
				    int nval, double* rhs)
{
  int i, j, node;
  double* xyz = param->xyz;
  int nnodes = this->GetNnodes();
  int ndofpn = this->GetNodeDOF();
  double d;
  //weights
  double* w = new double[nnodes];
  double wsum = 0.0;
  int nodeDOF;
  
  if(nval != ndofpn){
    std::cerr << "WARNING: # of values passed in does not match DOF per node!!" << std::endl;
    return;
  }

  //compute weights
  for(i = 0; i < nnodes; i++){
    d = 0.0;
    for(j = 0; j < 3; j++){
      node = this->nodes[i];
      d = (xyz[node*3 + j] - xyzPt[j])*(xyz[node*3 + j] - xyzPt[j]);
    }
    d = sqrt(d);
    w[i] = 1.0/d;
    wsum += w[i];
  }

  //Normalize weights
  for(i = 0; i < nnodes; i++){
    w[i] = w[i]/wsum;
  }

  //use inverse distance weighting to compute values at new the nodes
  for(j = 0; j < nnodes; j++){
    node = this->nodes[j];
    nodeDOF = param->nodeOffsetsDOF[node];
    for(i = 0; i < ndofpn; i++){
      rhs[nodeDOF + i] += values[i]/w[j];
    }
  }

  delete [] w;

  return;
}

}
