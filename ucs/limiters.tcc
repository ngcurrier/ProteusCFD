#include "eqnset_defines.h"
#include "driver.h"
#include "param.h"
#include "solutionSpace.h"
#include "solutionField.h"
#include "dataInfo.h"
#include <sstream>
#include <string>

template <class Type>
Limiter<Type>::Limiter(SolutionSpace<Type>* space, Type* q, Type* grad, Int neqn, Int stride, Int gradStride, std::string name)
{
  Int i;
  std::stringstream temposs;
  //set up internal data 
  this->nnodes = space->m->nnode;
  this->gnode = space->m->gnode;
  this->neqn = neqn;
  this->gradneqn = gradStride;
  this->nvars = stride;
  
  this->q = q;
  this->grad = grad;

  //setup data descriptor
  this->idata = new DataInfo(this->neqn, std::string("limiter_"+name));
  for(i = 0; i < neqn; i++){
    temposs.str("");
    temposs << i;
    this->idata->AddScalar(i, "limiter_"+name+"-"+temposs.str());
  }
  this->idata->Verify();

  //register the field on the solution space
  space->AddField(*this->idata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  l = space->GetField("limiter_"+name, FIELDS::STATE_NONE);
  
  //init limiters to one in case we never use them, otherwise we
  //could get junk values
  for(i = 0; i < (nnodes+gnode)*neqn; i++){
    l[i] = 1.0;
  }
  space->p->UpdateGeneralVectors(l, neqn);
}

template <class Type>
Limiter<Type>::~Limiter()
{
  delete this->idata;
}

template <class Type>
void Limiter<Type>::Compute(SolutionSpace<Type>* space)
{
  Int i, j;
  Param<Type>* param = space->param;

  //allocate qmin/qmax
  this->qmin = new Type[this->nnodes*this->neqn];
  this->qmax = new Type[this->nnodes*this->neqn];
  //initialize memory
  MemBlank(this->qmin, this->nnodes*this->neqn);
  MemBlank(this->qmax, this->nnodes*this->neqn);

  //initialize limiter
  for(i = 0; i < this->nnodes+this->gnode; i++){
    for(j = 0; j < this->neqn; j++){
      this->l[i*neqn + j] = 1.0;
    }
  }

  Kernel<Type> FindMinMax(Kernel_FindMinMax);
  Kernel<Type> BFindMinMax(Bkernel_FindMinMax);

  //call drivers to find local min/max
  DriverNoScatter(space, FindMinMax, this->neqn, (void*)this);
  BdriverNoScatter(space, BFindMinMax, this->neqn, (void*)this);
 
  //pick correct limiter
  if(param->limiter == 0){
    //shouldn't be here if limiters are turned off
  }
  else if(param->limiter == 1){
    Kernel<Type> BLimitKernel(Kernel_Barth);
    Kernel<Type> BBLimitKernel(Bkernel_Barth);
    //compute the limiter
    DriverNoScatter(space, BLimitKernel, this->neqn, (void*)this);
    BdriverNoScatter(space, BBLimitKernel, this->neqn, (void*)this);
  }
  else if(param->limiter == 2){
    Kernel<Type> VLimitKernel(Kernel_Venkat);
    Kernel<Type> VBLimitKernel(Bkernel_Venkat);
    //compute the limiter
    DriverNoScatter(space, VLimitKernel, this->neqn, (void*)this);
    BdriverNoScatter(space, VBLimitKernel, this->neqn, (void*)this);
  }
  else if(param->limiter == 3){
    Kernel<Type> VMLimitKernel(Kernel_VenkatMod);
    Kernel<Type> VMBLimitKernel(Bkernel_VenkatMod);
    //compute the limiter
    DriverNoScatter(space, VMLimitKernel, this->neqn, (void*)this);
    BdriverNoScatter(space, VMBLimitKernel, this->neqn, (void*)this);
  }
  else{
    //should never get here either
  }

  //check for negative pressures in eqnsets where that would be disastrous
  //only do this check if the q value we are looking at is actually the q
  //conservative solution variables, this limiter can also be run on turb
  //models, etc. and this doesn't make sense in that context
  if(param->limiter != 0){
    if(this->q == space->GetField("variableQ", FIELDS::STATE_NP1)){
      Kernel<Type> PressureClip(Kernel_PressureClip);
      DriverNoScatter(space, PressureClip, this->neqn, (void*)this);
    }
    //check all limiters for negative values and push them back to zero if found
    for(i = 0; i < this->nnodes; i++){
      for(j = 0; j < this->neqn; j++){
	if(real(this->l[i*neqn + j]) < 0.0){
	  l[i*neqn + j] = 0.0;
	}
      }
    }
  }
 
  //do a parallel sync of the limiter
  space->p->UpdateGeneralVectors(this->l, this->neqn);

  delete [] this->qmin;
  delete [] this->qmax;

  return;
}

template <class Type>
void Kernel_FindMinMax(KERNEL_ARGS)
{
  Int i;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;
  Type* qmaxL = limiter->qmax + neqn*left_cv;
  Type* qminL = limiter->qmin + neqn*left_cv;
  Type* qmaxR = limiter->qmax + neqn*right_cv;
  Type* qminR = limiter->qmin + neqn*right_cv;
  Type* q = limiter->q;
  Type* qL = q + nvars*left_cv;
  Type* qR = q + nvars*right_cv;
  
  for(i = 0; i < neqn; i++){
    qmaxL[i] = MAX(qmaxL[i], qR[i]);
    qminL[i] = MIN(qminL[i], qR[i]);
  }
  for(i = 0; i < neqn; i++){
    qmaxR[i] = MAX(qmaxR[i], qL[i]);
    qminR[i] = MIN(qminR[i], qL[i]);
  }
  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;

  return;
}

template <class Type>
void Bkernel_FindMinMax(B_KERNEL_ARGS)
{
  Int i;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Mesh<Type>* m = space->m;

  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;

  //only consider nodes for limiter if they are internal or parallel
  if(m->IsGhostNode(right_cv)){

    Type* qmaxL = limiter->qmax + neqn*left_cv;
    Type* qminL = limiter->qmin + neqn*left_cv;
    Type* q = limiter->q;
    Type* qR = q + nvars*right_cv;
    
    for(i = 0; i < neqn; i++){
      qmaxL[i] = MAX(qmaxL[i], qR[i]);
      qminL[i] = MIN(qminL[i], qR[i]);
    }
  }

  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;

  return;
}

template <class Type>
void Kernel_Barth(KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;
  Int gradneqn = limiter->gradneqn;
  Type* qmaxL = limiter->qmax + neqn*left_cv;
  Type* qminL = limiter->qmin + neqn*left_cv;
  Type* qmaxR = limiter->qmax + neqn*right_cv;
  Type* qminR = limiter->qmin + neqn*right_cv;
  Type* q = limiter->q;
  Type* qL = q + nvars*left_cv;
  Type* qR = q + nvars*right_cv;
  Type* xR,* xL;
  Type* gradL = limiter->grad + gradneqn*3*left_cv;
  Type* gradR = limiter->grad + gradneqn*3*right_cv;
  Type* limiterL = limiter->l + neqn*left_cv;
  Type* limiterR = limiter->l + neqn*right_cv;
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type temp;

  Type chi = param->chi;
  
  //compute extrapolations
  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;
  
  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);
  
  //TODO: generalize incase we want to use a non-midpoint edge CV
  dx[0] = 0.5*(xR[0] - xL[0]);
  dx[1] = 0.5*(xR[1] - xL[1]);
  dx[2] = 0.5*(xR[2] - xL[2]);
  
  for (j = 0; j < neqn; j++){
    QL[j] += 
      (0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
       (gradL[j*3 + 0]*dx[0] +
	gradL[j*3 + 1]*dx[1] +
	gradL[j*3 + 2]*dx[2]));
  }
  
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  
  for (j = 0; j < neqn; j++){
    QR[j] += 
      (0.5*chi*(qL[j] - qR[j]) + (1.0 - chi)*
       (gradR[j*3 + 0]*dx[0] +
	gradR[j*3 + 1]*dx[1] +
	gradR[j*3 + 2]*dx[2]));
  }
  
  //compute the Barth limiters now
  for(j = 0; j < neqn; j++){
    temp = 1.0;
    if(real(QL[j]) > real(qL[j])){
      temp = (qmaxL[j] - qL[j])/(QL[j] - qL[j]);
    }
    else if(real(QL[j]) < real(qL[j])){
      temp = (qminL[j] - qL[j])/(QL[j] - qL[j]);
    }
    
    //watch out for negative numbers
    temp = MAX(0.0, temp);
    //compute the maximum allowable limiter
    temp = MIN(1.0, temp);
    
    //pick the smallest phi required
    limiterL[j] = MIN(limiterL[j], temp);
  }
  
  for(j = 0; j < neqn; j++){
    temp = 1.0;
    if(real(QR[j]) > real(qR[j])){
      temp = (qmaxR[j] - qR[j])/(QR[j] - qR[j]);
    }
    else if(real(QR[j]) < real(qR[j])){
      temp = (qminR[j] - qR[j])/(QR[j] - qR[j]);
    }
    
    //watch out for negative numbers
    temp = MAX(0.0, temp);
    //compute the maximum allowable limiter
    temp = MIN(1.0, temp);
    
    //pick the smallest phi required
    limiterR[j] = MIN(limiterR[j], temp);
  }
  

  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;
  
  return;
}


template <class Type>
void Bkernel_Barth(B_KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Int neqn = limiter->neqn;
    Int nvars = limiter->nvars;
    Int gradneqn = limiter->gradneqn;
    Type* qmaxL = limiter->qmax + neqn*left_cv;
    Type* qminL = limiter->qmin + neqn*left_cv;
    Type* q = limiter->q;
    Type* qL = q + nvars*left_cv;
    Type* qR = q + nvars*right_cv;
    Type* xR,* xL;
    Type* gradL = limiter->grad + gradneqn*3*left_cv;
    Type* limiterL = limiter->l + neqn*left_cv;
    Type* QL = (Type*)alloca(sizeof(Type)*nvars);
    Type* QR = (Type*)alloca(sizeof(Type)*nvars);
    Type* dx = (Type*)alloca(sizeof(Type)*3);
    Type temp;
    
    Type chi = param->chi;
    
    //compute extrapolations
    xL = m->cg + 3*left_cv;
    xR = m->cg + 3*right_cv;
    
    memcpy(QL, qL, sizeof(Type)*nvars);
    memcpy(QR, qR, sizeof(Type)*nvars);
    
    //TODO: generalize incase we want to use a non-midpoint edge CV
    dx[0] = 0.5*(xR[0] - xL[0]);
    dx[1] = 0.5*(xR[1] - xL[1]);
    dx[2] = 0.5*(xR[2] - xL[2]);
    
    for (j = 0; j < neqn; j++){
      QL[j] += 
	(0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
	 (gradL[j*3 + 0]*dx[0] +
	  gradL[j*3 + 1]*dx[1] +
	  gradL[j*3 + 2]*dx[2]));
    }
    
    //compute the Barth limiters now
    for(j = 0; j < neqn; j++){
      temp = 1.0;
      if(real(QL[j]) > real(qL[j])){
	temp = (qmaxL[j] - qL[j])/(QL[j] - qL[j]);
      }
      else if(real(QL[j]) < real(qL[j])){
	temp = (qminL[j] - qL[j])/(QL[j] - qL[j]);
      }
      
      //watch out for negative numbers
      temp = MAX(0.0, temp);
      //compute the maximum allowable limiter
      temp = MIN(1.0, temp);
      
      //pick the smallest phi required
      limiterL[j] = MIN(limiterL[j], temp);
    }
  }
  
  //no scatter
  *size = 0;
  *ptrL = NULL;
  
  return;
}

template <class Type>
void Kernel_Venkat(KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Mesh<Type>* m = space->m;
  Param<Type>* param = space->param;
  
  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;
  Int gradneqn = limiter->gradneqn;
  Type* qmaxL = limiter->qmax + neqn*left_cv;
  Type* qminL = limiter->qmin + neqn*left_cv;
  Type* qmaxR = limiter->qmax + neqn*right_cv;
  Type* qminR = limiter->qmin + neqn*right_cv;
  Type* q = limiter->q;
  Type* qL = q + nvars*left_cv;
  Type* qR = q + nvars*right_cv;
  Type* xR,* xL;
  Type* gradL = limiter->grad + gradneqn*3*left_cv;
  Type* gradR = limiter->grad + gradneqn*3*right_cv;
  Type* limiterL = limiter->l + neqn*left_cv;
  Type* limiterR = limiter->l + neqn*right_cv;
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type temp;

  Type chi = param->chi;
  
  //compute extrapolations
  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;
  
  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);
  
  //TODO: generalize incase we want to use a non-midpoint edge CV
  dx[0] = 0.5*(xR[0] - xL[0]);
  dx[1] = 0.5*(xR[1] - xL[1]);
  dx[2] = 0.5*(xR[2] - xL[2]);
  
  for (j = 0; j < neqn; j++){
    QL[j] += 
      (0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
       (gradL[j*3 + 0]*dx[0] +
	gradL[j*3 + 1]*dx[1] +
	gradL[j*3 + 2]*dx[2]));
  }
  
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  
  for (j = 0; j < neqn; j++){
    QR[j] += 
      (0.5*chi*(qL[j] - qR[j]) + (1.0 - chi)*
       (gradR[j*3 + 0]*dx[0] +
	gradR[j*3 + 1]*dx[1] +
	gradR[j*3 + 2]*dx[2]));
  }
  
  //compute the Venkat limiters now
  for(j = 0; j < neqn; j++){
    temp = 1.0;
    if(real(QL[j]) > real(qL[j])){
      temp = (qmaxL[j] - qL[j])/(QL[j] - qL[j]);
    }
    else if(real(QL[j]) < real(qL[j])){
      temp = (qminL[j] - qL[j])/(QL[j] - qL[j]);
    }
    
    temp = (temp*temp + 2.0*temp)/(temp*temp + temp + 2.0);
    limiterL[j] = MIN(limiterL[j], temp);
  }
  
  for(j = 0; j < neqn; j++){
    temp = 1.0;
    if(real(QR[j]) > real(qR[j])){
      temp = (qmaxR[j] - qR[j])/(QR[j] - qR[j]);
    }
    else if(real(QR[j]) < real(qR[j])){
      temp = (qminR[j] - qR[j])/(QR[j] - qR[j]);
    }
    
    temp = (temp*temp + 2.0*temp)/(temp*temp + temp + 2.0);
    limiterR[j] = MIN(limiterR[j], temp);
  }
  

  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;
  
  return;
}


template <class Type>
void Bkernel_Venkat(B_KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Int neqn = limiter->neqn;
    Int nvars = limiter->nvars;
    Int gradneqn = limiter->gradneqn;
    Type* qmaxL = limiter->qmax + neqn*left_cv;
    Type* qminL = limiter->qmin + neqn*left_cv;
    Type* q = limiter->q;
    Type* qL = q + nvars*left_cv;
    Type* qR = q + nvars*right_cv;
    Type* xR,* xL;
    Type* gradL = limiter->grad + gradneqn*3*left_cv;
    Type* limiterL = limiter->l + neqn*left_cv;
    Type* QL = (Type*)alloca(sizeof(Type)*nvars);
    Type* QR = (Type*)alloca(sizeof(Type)*nvars);
    Type* dx = (Type*)alloca(sizeof(Type)*3);
    Type temp;
    
    Type chi = param->chi;
    
    //compute extrapolations
    xL = m->cg + 3*left_cv;
    xR = m->cg + 3*right_cv;
    
    memcpy(QL, qL, sizeof(Type)*nvars);
    memcpy(QR, qR, sizeof(Type)*nvars);
    
    //TODO: generalize incase we want to use a non-midpoint edge CV
    dx[0] = 0.5*(xR[0] - xL[0]);
    dx[1] = 0.5*(xR[1] - xL[1]);
    dx[2] = 0.5*(xR[2] - xL[2]);
    
    for (j = 0; j < neqn; j++){
      QL[j] += 
	(0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
	 (gradL[j*3 + 0]*dx[0] +
	  gradL[j*3 + 1]*dx[1] +
	  gradL[j*3 + 2]*dx[2]));
    }
    
    //compute the Venkat limiters now
    for(j = 0; j < neqn; j++){
      temp = 1.0;
      if(real(QL[j]) > real(qL[j])){
	temp = (qmaxL[j] - qL[j])/(QL[j] - qL[j]);
      }
      else if(real(QL[j]) < real(qL[j])){
	temp = (qminL[j] - qL[j])/(QL[j] - qL[j]);
      }
      
      temp = (temp*temp + 2.0*temp)/(temp*temp + temp + 2.0);
      limiterL[j] = MIN(limiterL[j], temp);
    }
  }
  
  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;
  
  return;
}

template <class Type>
void Kernel_VenkatMod(KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;
  Int gradneqn = limiter->gradneqn;
  Type* qmaxL = limiter->qmax + neqn*left_cv;
  Type* qminL = limiter->qmin + neqn*left_cv;
  Type* qmaxR = limiter->qmax + neqn*right_cv;
  Type* qminR = limiter->qmin + neqn*right_cv;
  Type* q = limiter->q;
  Type* qL = q + nvars*left_cv;
  Type* qR = q + nvars*right_cv;
  Type* xR,* xL;
  Type* gradL = limiter->grad + gradneqn*3*left_cv;
  Type* gradR = limiter->grad + gradneqn*3*right_cv;
  Type* limiterL = limiter->l + neqn*left_cv;
  Type* limiterR = limiter->l + neqn*right_cv;
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type temp;

  Type chi = param->chi;

  //Pi
  Type Pi = 3.141592653589793;

  //tunable parameters
  Type K = 1.0;
  Type K3 = K*K*K;

  //compute extrapolations
  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;
  
  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);

  //compute a characteristic length of the control volume
  //use a sphere diameter with the same volume 4/3*pi*r^3
  Type Ll3 = m->vol[left_cv];  
  Type Lr3 = m->vol[right_cv];
  Type ep2;
  Type DP, DM;

  Ll3 = 6.0*Pi*Ll3;
  Lr3 = 6.0*Pi*Lr3;

  //TODO: generalize incase we want to use a non-midpoint edge CV
  dx[0] = 0.5*(xR[0] - xL[0]);
  dx[1] = 0.5*(xR[1] - xL[1]);
  dx[2] = 0.5*(xR[2] - xL[2]);
  
  for (j = 0; j < neqn; j++){
    QL[j] += 
      (0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
       (gradL[j*3 + 0]*dx[0] +
	gradL[j*3 + 1]*dx[1] +
	gradL[j*3 + 2]*dx[2]));
  }
  
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  
  for (j = 0; j < neqn; j++){
    QR[j] += 
      (0.5*chi*(qL[j] - qR[j]) + (1.0 - chi)*
       (gradR[j*3 + 0]*dx[0] +
	gradR[j*3 + 1]*dx[1] +
	gradR[j*3 + 2]*dx[2]));
  }
  
  //compute the Venkat limiters now
  for(j = 0; j < neqn; j++){
    //compute delta+ and delta- for the limiter
    DM = QL[j] - qL[j];
    if(real(QL[j]) > real(qL[j])){
      DP = qmaxL[j] - qL[j];
    }
    else{
      DP = qminL[j] - qL[j];
    }
    
    //compute eps^2
    ep2 = Ll3*K3;
    temp = (DP*DP + ep2 + 2.0*DM*DP)/(DP*DP + 2.0*DM*DM + DM*DP + ep2);
    limiterL[j] = MIN(limiterL[j], temp);
  }
  
  for(j = 0; j < neqn; j++){
    //compute delta+ and delta- for the limiter
    DM = QR[j] - qR[j];
    if(real(QR[j]) > real(qR[j])){
      DP = qmaxR[j] - qR[j];
    }
    else{
      DP = qminR[j] - qR[j];
    }

    //compute eps^2
    ep2 = Lr3*K3;
    temp = (DP*DP + ep2 + 2.0*DM*DP)/(DP*DP + 2.0*DM*DM + DM*DP + ep2);
    limiterR[j] = MIN(limiterR[j], temp);
  }
  

  //no scatter
  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;
  
  return;
}


template <class Type>
void Bkernel_VenkatMod(B_KERNEL_ARGS)
{
  Int j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Int neqn = limiter->neqn;
    Int nvars = limiter->nvars;
    Int gradneqn = limiter->gradneqn;
    Type* qmaxL = limiter->qmax + neqn*left_cv;
    Type* qminL = limiter->qmin + neqn*left_cv;
    Type* q = limiter->q;
    Type* qL = q + nvars*left_cv;
    Type* qR = q + nvars*right_cv;
    Type* xR,* xL;
    Type* gradL = limiter->grad + gradneqn*3*left_cv;
    Type* limiterL = limiter->l + neqn*left_cv;
    Type* QL = (Type*)alloca(sizeof(Type)*nvars);
    Type* QR = (Type*)alloca(sizeof(Type)*nvars);
    Type* dx = (Type*)alloca(sizeof(Type)*3);
    Type temp;

    Type chi = param->chi;
    
    //Pi
    Type Pi = 3.141592653589793;

    //tunable parameters
    Type K = 1.0;
    Type K3 = K*K*K;

    //compute a characteristic length of the control volume
    //use a sphere diameter with the same volume 4/3*pi*r^3
    Type Ll3 = m->vol[left_cv];  
    Ll3 = 6.0*Pi*Ll3;
    
    Type DM, DP, ep2;
    
    //compute extrapolations
    xL = m->cg + 3*left_cv;
    xR = m->cg + 3*right_cv;
    
    memcpy(QL, qL, sizeof(Type)*nvars);
    memcpy(QR, qR, sizeof(Type)*nvars);
    
    //TODO: generalize incase we want to use a non-midpoint edge CV
    dx[0] = 0.5*(xR[0] - xL[0]);
    dx[1] = 0.5*(xR[1] - xL[1]);
    dx[2] = 0.5*(xR[2] - xL[2]);
    
    for (j = 0; j < neqn; j++){
      QL[j] += 
	(0.5*chi*(qR[j] - qL[j]) + (1.0 - chi)*
	 (gradL[j*3 + 0]*dx[0] +
	  gradL[j*3 + 1]*dx[1] +
	  gradL[j*3 + 2]*dx[2]));
    }
    
    //compute the Venkat limiters now
    for(j = 0; j < neqn; j++){
      //compute delta+ and delta- for the limiter
      DM = QL[j] - qL[j];
      if(real(QL[j]) > real(qL[j])){
	DP = qmaxL[j] - QL[j];
      }
      else{
	DP = qminL[j] - QL[j];
      }
      
      //compute eps^2
      ep2 = Ll3*K3;
      temp = (DP*DP + ep2 + 2.0*DM*DP)/(DP*DP + 2.0*DM*DM + DM*DP + ep2);
      limiterL[j] = MIN(limiterL[j], temp);
    }
  }
  
  //no scatter
  *size = 0;
  *ptrL = NULL;
  
  return;
}

template <class Type>
void Kernel_PressureClip(KERNEL_ARGS)
{
  Int i, j;
  Limiter<Type>* limiter = (Limiter<Type>*) custom;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;

  Int neqn = limiter->neqn;
  Int nvars = limiter->nvars;
  Int gradneqn = limiter->gradneqn;
  Type* q = limiter->q;
  Type* qL = q + nvars*left_cv;
  Type* qR = q + nvars*right_cv;
  Type* xR,* xL;
  Type* gradL = limiter->grad + gradneqn*3*left_cv;
  Type* gradR = limiter->grad + gradneqn*3*right_cv;
  Type* limiterL = limiter->l + neqn*left_cv;
  Type* limiterR = limiter->l + neqn*right_cv;
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* Qroe = (Type*)alloca(sizeof(Type)*nvars);
  Type* dQ = (Type*)alloca(sizeof(Type)*neqn);
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type PL, PR;
  Type Pmin = 1.0e-10;
  Type Emin = 1.0e-10;

  Type chi = param->chi;
  Type gamma = param->gamma;

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;

  //TODO: generalize incase we want to use a non-midpoint edge CV
  dx[0] = 0.5*(xR[0] - xL[0]);
  dx[1] = 0.5*(xR[1] - xL[1]);
  dx[2] = 0.5*(xR[2] - xL[2]);
  
  for (j = 0; j < neqn; j++){
    dQ[j] = qR[j] - qL[j];
  }
  eqnset->ExtrapolateVariables(QL, qL, dQ, gradL, dx, limiterL);

  //check for sanity of extrapolation
  if(eqnset->BadExtrapolation(QL)){
    for(i = 0; i < neqn; i++){
      limiterL[i] = 0.0;
    }
  }

  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  
  for (j = 0; j < neqn; j++){
    dQ[j] = -dQ[j];
  }
  eqnset->ExtrapolateVariables(QR, qR, dQ, gradR, dx, limiterR);

  //check for sanity of extrapolation
  if(eqnset->BadExtrapolation(QR)){
    for(i = 0; i < neqn; i++){
      limiterR[i] = 0.0;
    }
  }

  //also check the roe variables call, if pressure turns out to be negative judo chop
  //both sides, since not all eqnsets implement RoeVariables, it will return true only
  //if the output will be meaningful
  if(eqnset->RoeVariables(QL, QR, gamma, Qroe)){
    eqnset->ComputeAuxiliaryVariables(Qroe);
    if(eqnset->BadExtrapolation(Qroe)){
      for(i = 0; i < neqn; i++){
	limiterL[i] = limiterR[i] = 0.0;
      }
    }    
  }

  return;
}
