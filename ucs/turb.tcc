#include "mesh.h"
#include "driver.h"
#include "bc.h"
#include "gradient.h"
#include "jacobian.h"
#include "solutionSpace.h"
#include "exceptions.h"

template <class Type>
TurbulenceModel<Type>::TurbulenceModel(SolutionSpace<Type>* space):
  neqn(0), space(space), tvar(NULL), tvarinf(NULL), tgrad(NULL), idata(NULL), limiter(NULL)
{
  //add turbulent viscosity field
  space->AddField("mut");
}

template <class Type>
TurbulenceModel<Type>::~TurbulenceModel()
{
  delete [] tvarinf;
  delete limiter;
  delete idata;

  return;
}

template <class Type>
void TurbulenceModel<Type>::SetTinf()
{
  //do nothing by default
  return;
}

template <class Type>
void TurbulenceModel<Type>::Initialize()
{
  Abort << "WARNING: In turbulence model Initialize() not implemented!!!";
  return;
}

template <class Type>
void TurbulenceModel<Type>::UpdateBCs()
{
  Kernel<Type> BC(BC_Turb_Kernel);

  //call driver to loop over and set boundary conditions
  BdriverNoScatter(space, BC, this->neqn, this);

  return;
}

template <class Type>
void TurbulenceModel<Type>::BC_Kernel(B_KERNEL_ARGS)
{
  return;
}

template <class Type>
void TurbulenceModel<Type>::BC_Jac_Kernel(B_KERNEL_ARGS)
{
  return;
}

template <class Type>
void TurbulenceModel<Type>::Source(Type nu, Type d, Type* vgrad, Type* tvars, Type vol, 
				   Type* res, Type* jac)
{
  Abort << "WARNING: In turbulence model Source() not implemented!!!";
  return;
}

template <class Type>
void TurbulenceModel<Type>::DiffusiveDriver()
{
  Kernel<Type> Diffusive_Flux(Kernel_Diffusive);
  Kernel<Type> Bdiffusive_Flux(Bkernel_Diffusive);

  //call drivers to compute the diffusive part of the flux
  DriverNoScatter(space, Diffusive_Flux, neqn, this);
  BdriverNoScatter(space, Bdiffusive_Flux, neqn, this);

  return;
}

template <class Type>
void TurbulenceModel<Type>::Diffusive(Type nu, Type* tgrad, Type* tvarsL, Type* tvarsR, 
				      Type* avec, Type dgrad, Type* resL, Type* resR, 
				      Type* jacL, Type* jacR)
{
  Abort << "WARNING: In turbulence model Diffusive() not implemented!!!";
  return;
}

template <class Type>
void TurbulenceModel<Type>::Convective()
{
  Kernel<Type> Convective_Flux(Kernel_Convective);
  Kernel<Type> Bconvective_Flux(Bkernel_Convective);

  //call drivers to compute the convective part of the flux
  Driver(space, Convective_Flux, neqn, this);
  Bdriver(space, Bconvective_Flux, neqn, this);
    
  return;
}

template <class Type>
Type TurbulenceModel<Type>::ComputeEddyViscosity(Type rho, Type nu, Int node)
{
  //std::cerr << "WARNING: In turbulence model ComputeEddyViscosity() not implemented!!!" << std::endl;
  return 0.0;
}

template <class Type>
void TurbulenceModel<Type>::TemporalResidual(Type timestep, Int iter, Int torder)
{
  Int offset;
  Type dtm1, dt;
  Type cnp1, cnm1;
  Type* dq = new Type[neqn];
  Type* dqm1 = new Type[neqn];
  Mesh<Type>* m = space->m;
  Int nnode = m->GetNumNodes();

  if(iter > 1 && torder == 2){
    //coefficients for BDF2 are phi_n+1 = 1.5, phi_n = -2.0, phi_n-1 = 0.5 
    cnp1 = 1.5;
    cnm1 = -0.5;
    //cn = -2.0;
  }
  else{
    cnp1 = 1.0;
    cnm1 = 0.0;
    //cn = -1.0;
  }
  for(Int i = 0; i < nnode; i++){
    offset = i*neqn;
    dt = cnp1*m->vol[i]/timestep;           //use v_n+1
    if(space->eqnset->param->gcl && iter > 2){
      dtm1 = cnm1*m->vololdm1[i]/timestep;  //use v_n-1
    }
    else{
      dtm1 = cnm1*m->vol[i]/timestep;       //use v_n+1
    }
    Type* res = &crs.b[i*neqn];
    for(Int j = 0; j < neqn; j++){
      dq[j] = tvar[offset + j] - tvarold[offset + j];
      dqm1[j] = tvarold[offset + j] - tvaroldm1[offset + j];
    }
    for(Int j = 0; j < neqn; j++){
      res[j] -= dt*dq[j];
      res[j] -= dtm1*dqm1[j];
    }
  }
  delete [] dq;
  delete [] dqm1;
}

template <class Type>
void TurbulenceModel<Type>::ContributeTemporalTerms(Type vol, Type cnp1, Type dt, Type dtau, Type* A)
{
  for(Int i = 0; i < neqn; i++){
    if(space->param->useLocalTimeStepping && (real(space->param->dt) > 0.0)){
      A[i*neqn + i] += cnp1*vol/dt + vol/dtau;
    }
    else{
      A[i*neqn + i] += cnp1*vol/dtau;
    }
  }
}

template <class Type>
void TurbulenceModel<Type>::Compute()
{
  Int i, j;
  Mesh<Type>* m = space->m;
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  EqnSet<Type>* eqnset = space->eqnset;
  Type* tres = (Type*)alloca(sizeof(Type)*neqn);
  Type* tjac = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Int qneqn = eqnset->neqn;
  Int qnvars = qneqn + eqnset->nauxvars;
  Param<Type>* param = space->param;
  Type* dt = space->GetFieldData("timestep", FIELDS::STATE_NONE);
  Type* mut = space->GetFieldData("mut", FIELDS::STATE_NONE);

  //blank the crs system
  crs.BlankSystem();
  
  //compute new bcs
  UpdateBCs();
  //perform parallel update
  space->p->UpdateGeneralVectors(tvar, neqn);

  //compute gradients for the turbulent variables
  Int weighted = true;
  Gradient<Type> gradt(this->neqn, this->neqn, NULL, tvar, space, param->gradType, 
		       tgrad, weighted);
  gradt.Compute();
  
  //if we are doing higher order compute limiter
  if((param->turbModelSorder > 1) && (space->iter > param->nFirstOrderSteps)){
    limiter->Compute(space);
  }

  //compute convective terms
  Convective();

  //compute diffusive terms
  DiffusiveDriver();

  std::vector<Int> gradientLoc;
  Int nterms = eqnset->GetGradientsLocation(gradientLoc);
  Int vloc = eqnset->GetVelocityGradLocation()*3;
  
  //compute source terms
  Type* dist = space->GetFieldData("wallDistance", FIELDS::STATE_NONE);
  for(i = 0; i < nnode; i++){
    Type* q = &space->q[i*qnvars];
    Type vol = m->vol[i];
    Type d = dist[i];
    Type rho = eqnset->GetDensity(q);
    Type mu = eqnset->ComputeViscosity(q);
    Type nu = mu/rho;
    Type* res = &crs.b[i*neqn];
    Type* jacd = crs.A->GetPointer(i, i);
    Type* qGrad = &space->qgrad[i*nterms*3];
    Type* vgrad = &qGrad[vloc];
    Source(nu, d, vgrad, &tvar[i*neqn], vol, tres, tjac);
    for(j = 0; j < neqn; j++){
      res[j] += tres[j];
    }
    for(j = 0; j < neqn*neqn; j++){
      jacd[j] += tjac[j];
    }
  }

  //compute temporal derivative terms
  if(param->torder){
    TemporalResidual(param->dt, space->iter, param->torder);
  }

  //sum to the diagonal
  Kernel<Type> Diag_NumJac(Kernel_Diag_NumJac);
  Driver(space, Diag_NumJac, neqn*neqn, (void*)&crs);

  //contribute to the diagonal jacobian
  Type cnp1 = 1.0;
  if(space->iter > 1 && param->torder == 2){
    cnp1 = 1.5;
  }
  for(Int i = 0; i < nnode; i++){
    ContributeTemporalTerms(m->vol[i], cnp1, param->dt, dt[i], crs.A->GetPointer(i,i));
  }

  //modify the CRS system if necessary
  Kernel<Type> Jac_Modify(BC_Turb_Jac_Kernel);
  BdriverNoScatter(space, Jac_Modify, neqn*neqn, this);

  Type resid = VecL2Norm(crs.b, nnode*neqn);
  Type residGlobal = ParallelL2Norm(space->p, crs.b, nnode*neqn);

  std::cout << "||Turb res||: " << resid << " ";

  if(space->p->GetRank() == 0){
    space->residOutFile << "||Turb-res||: " << residGlobal << " ";
  }

  if(param->nSgs > 0){
    crs.A->PrepareSGS();
    Int nsgs = MAX(param->nSgs, -param->nSgs);
    //report the true residual here
    Type solveres = crs.SGS(nsgs, NULL, NULL, NULL, 0);
    std::cout << "||Turb-Dq||: " << solveres << " ";
    if(std::isnan(real(solveres)) || std::isinf(real(solveres))){
      Abort << "Turbulence solution residual divergent!! Infinite residual -- I'm OUT!";
    }
  }
  else{
    for(i = 0; i < nnode; i++){
      for(j = 0; j < neqn; j++){
	crs.x[i*neqn + j] = crs.b[i*neqn + j]*dt[i]/m->vol[i];
      }
    }
  }
  

  //now apply the solution change
  Type relax = 0.5;
  for(i = 0; i < nnode; i++){
    for(j = 0; j < neqn; j++){
      //tvar[i*neqn + j] += MIN(tvar[i*neqn + j]*relax, crs.x[i*neqn + j]);
      tvar[i*neqn + j] += crs.x[i*neqn + j];
      if(real(tvar[i*neqn + j]) < 0.0){
	std::cerr << "clipping turb var dtq: " << crs.x[i] << " tq: " << tvar[i] << std::endl;
	std::cerr << "\t At coordinates: " << m->xyz[i*3 + 0] << " " << m->xyz[i*3 + 1] 
		  << " " << m->xyz[i*3 + 2] << std::endl;
	tvar[i*neqn + j] = 0.0;
      }
    }
  }
  //perform parallel update
  space->p->UpdateGeneralVectors(tvar, neqn);

  //now that we have computed the turbulence variables
  //compute the eddy viscosity which is what we actually need
  for(i = 0; i < nnode+gnode; i++){
    Type* q = &space->q[i*qnvars];
    Type rho = eqnset->GetDensity(q);
    Type T = eqnset->GetTemperature(q);
    Type mu = eqnset->ComputeViscosity(q);
    Type nu = mu/rho;
    mut[i] = ComputeEddyViscosity(rho, nu, i);
  }

  return;
}


template <class Type>
void Kernel_Convective(KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;

  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  Int neqn = turb->neqn;
  Int qneqn = eqnset->neqn;
  Int qnvars = qneqn + eqnset->nauxvars;
  Type chi = param->chi;

  Type* qL = &space->q[left_cv*qnvars];
  Type* qR = &space->q[right_cv*qnvars];
  
  Type* xL = &m->cg[left_cv*3];
  Type* xR = &m->cg[right_cv*3];
  Type* dx = (Type*)alloca(sizeof(Type)*3);

  Type* q = (Type*)alloca(sizeof(Type)*qneqn);
  Type* TL = (Type*)alloca(sizeof(Type)*neqn);
  Type* TR = (Type*)alloca(sizeof(Type)*neqn);

  Type* limiterL = &turb->limiter->l[left_cv*neqn];
  Type* limiterR = &turb->limiter->l[right_cv*neqn];
  
  for(i = 0; i < qneqn; i++){
    q[i] = 0.5*(qL[i] + qR[i]);
  }

  Type* tL = &turb->tvar[left_cv*neqn];
  Type* tR = &turb->tvar[right_cv*neqn];
  
  //copy turb vars for h.o. extrapolation
  memcpy(TL, tL, sizeof(Type)*neqn);
  memcpy(TR, tR, sizeof(Type)*neqn);

  //use strict upwinding
  Type theta = eqnset->GetTheta(q, avec, vdotn);
  Type area = avec[3];

  Type* jacR = turb->crs.A->GetPointer(left_cv, right_cv);
  Type* jacL = turb->crs.A->GetPointer(right_cv, left_cv);

  if(param->turbModelSorder > 1 && (space->iter > param->nFirstOrderSteps)){
    Type* tgradL = &turb->tgrad[left_cv*neqn*3];
    Type* tgradR = &turb->tgrad[right_cv*neqn*3];

    //TODO: generalize incase we want to use a non-midpoint edge CV
    dx[0] = 0.5*(xR[0] - xL[0]);
    dx[1] = 0.5*(xR[1] - xL[1]);
    dx[2] = 0.5*(xR[2] - xL[2]);
    for (j = 0; j < neqn; j++){
      TL[j] += 
	(0.5*chi*(tR[j] - tL[j]) + (1.0 - chi)*
	 (tgradL[j*3 + 0]*dx[0] +
	  tgradL[j*3 + 1]*dx[1] +
	  tgradL[j*3 + 2]*dx[2]))*limiterL[j];
    }
    
    dx[0] = -dx[0];
    dx[1] = -dx[1];
    dx[2] = -dx[2];
    
    for (j = 0; j < neqn; j++){
      TR[j] += 
	(0.5*chi*(tL[j] - tR[j]) + (1.0 - chi)*
	 (tgradR[j*3 + 0]*dx[0] +
	  tgradR[j*3 + 1]*dx[1] +
	  tgradR[j*3 + 2]*dx[2]))*limiterR[j];
    }
  } 
 //flux from left to right
  if(real(theta) > 0.0){
    for(i = 0; i < neqn; i++){
      //compute the derivative
      tempR[i] = theta*area;
      //contribute to the left jacobian
      jacL[i*neqn + i] -= tempR[i];
      //finish computing the convective residual
      tempR[i] *= TL[i];
    }
  }
  //flux from right to left
  else{
    for(i = 0; i < neqn; i++){
      //compute the derivative
      tempR[i] = theta*area;
      //contribute to the right jacobian
      jacR[i*neqn + i] += tempR[i];
      //finish compute the convective residual
      tempR[i] *= TR[i];
    }
  }

  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempR[i];
  }

  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &turb->crs.b[left_cv*neqn];
  *ptrR = &turb->crs.b[right_cv*neqn];
  
  return;
}

template <class Type>
void Bkernel_Convective(B_KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;
  Mesh<Type>* m = space->m;
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  Int neqn = turb->neqn;
  Int qneqn = eqnset->neqn;
  Int qnvars = qneqn + eqnset->nauxvars;

  Type chi = param->chi;

  Type* qL = &space->q[left_cv*qnvars];
  Type* qR = &space->q[right_cv*qnvars];
  
  Type* xL = &m->cg[left_cv*3];
  Type* xR = &m->cg[right_cv*3];
  Type* dx = (Type*)alloca(sizeof(Type)*3);

  Type* q = (Type*)alloca(sizeof(Type)*qneqn);
  Type* TL = (Type*)alloca(sizeof(Type)*neqn);
  Type* TR = (Type*)alloca(sizeof(Type)*neqn);

  for(i = 0; i < qneqn; i++){
    q[i] = 0.5*(qL[i] + qR[i]);
  }

  Type* tL = &turb->tvar[left_cv*neqn];
  Type* tR = &turb->tvar[right_cv*neqn];

  //copy turb vars for h.o. extrapolation
  memcpy(TL, tL, sizeof(Type)*neqn);
  memcpy(TR, tR, sizeof(Type)*neqn);

  Type* limiterL = &turb->limiter->l[left_cv*neqn];
  
  //use strict upwinding
  Type theta = eqnset->GetTheta(q, avec, vdotn);
  Type area = avec[3];

  if(param->turbModelSorder > 1 && (space->iter > param->nFirstOrderSteps)){
    //only do full extrapolation if boundary node is ghost(parallel)
    if(m->IsGhostNode(right_cv)){
      Type* tgradL = &turb->tgrad[left_cv*neqn*3];
      
      //TODO: generalize incase we want to use a non-midpoint edge CV
      dx[0] = 0.5*(xR[0] - xL[0]);
      dx[1] = 0.5*(xR[1] - xL[1]);
      dx[2] = 0.5*(xR[2] - xL[2]);
      for (j = 0; j < neqn; j++){
	TL[j] += 
	  (0.5*chi*(tR[j] - tL[j]) + (1.0 - chi)*
	   (tgradL[j*3 + 0]*dx[0] +
	    tgradL[j*3 + 1]*dx[1] +
	    tgradL[j*3 + 2]*dx[2]))*limiterL[j];
      }
      dx[0] = -dx[0];
      dx[1] = -dx[1];
      dx[2] = -dx[2];
      
      Type* tgradR = &turb->tgrad[right_cv*neqn*3];
      Type* limiterR = &turb->limiter->l[right_cv*neqn];
      for (j = 0; j < neqn; j++){
	TR[j] += 
	  (0.5*chi*(tL[j] - tR[j]) + (1.0 - chi)*
	   (tgradR[j*3 + 0]*dx[0] +
	    tgradR[j*3 + 1]*dx[1] +
	    tgradR[j*3 + 2]*dx[2]))*limiterR[j];
      }
    }
  }
  
  Type* diagL = turb->crs.A->GetPointer(left_cv, left_cv);
  //flux from left to right
  if(real(theta) > 0.0){
    for(i = 0; i < neqn; i++){
      //compute the derivative
      tempR[i] = theta*area;
      //flux leaving
      diagL[i*neqn + i] += tempR[i];
      //finish computing the convective residual
      tempR[i] *= TL[i];
    }
  }
  //flux from right to left
  else{
    for(i = 0; i < neqn; i++){
      //compute the derivative
      tempR[i] = theta*area;
      //this is if we have a parallel node on the boundary edge, i.e. in the volume
      if(m->IsGhostNode(right_cv)){
	Type* jacR = turb->crs.A->GetPointer(left_cv, right_cv);
	jacR[i*neqn + i] += tempR[i];
      }
      tempR[i] *= TR[i];
    }
  }
  
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempR[i];
  }

  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &turb->crs.b[left_cv*neqn];
  *ptrR = NULL;

  return;
}

template <class Type>
void Kernel_Diffusive(KERNEL_ARGS)
{  
  Int i, j;
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  Int neqn = turb->neqn;

  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;

  Type* tvarsL = &turb->tvar[left_cv*neqn];
  Type* tvarsR = &turb->tvar[right_cv*neqn];

  Type* tresL = (Type*)alloca(sizeof(Type)*neqn);
  Type* tresR = (Type*)alloca(sizeof(Type)*neqn);
  Type* tjacL = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* tjacR = (Type*)alloca(sizeof(Type)*neqn*neqn);

  Type dgrad;
  Type ds2;
  Type* de = (Type*)alloca(sizeof(Type)*3);

  ds2 = 0.0;
  for(i = 0; i < 3; i++){
    de[i] = m->cg[right_cv*3 + i] - m->cg[left_cv*3 + i];
    ds2 += de[i]*de[i];
  }
  
  Type dx, dy, dz, d;
  dx = de[0]*avec[0];
  dy = de[1]*avec[1];
  dz = de[2]*avec[2];

  d = dx + dy + dz;
  dgrad = d/ds2;

  Int qneqn = eqnset->neqn;
  Int qnvars = qneqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*qnvars];
  Type* qR = &space->q[right_cv*qnvars];
  Type* qavg = (Type*)alloca(sizeof(Type)*qnvars);

  for(i = 0; i < qneqn; i++){
    qavg[i] = 0.5*(qL[i] + qR[i]);
  }
  eqnset->ComputeAuxiliaryVariables(qavg);
  Type rho = eqnset->GetDensity(qavg);
  Type mu = eqnset->ComputeViscosity(qavg);
  Type nu = mu/rho;

  Type* tgrad = (Type*)alloca(sizeof(Type)*neqn*3);
  for(i = 0; i < neqn*3; i++){
    tgrad[i] = 0.5*(turb->tgrad[left_cv*neqn*3 + i] + turb->tgrad[right_cv*neqn*3 + i]);
  }
  
  //do the h.o. extrapolation (directional derivative)
  Type qdots, dq;
  for(j = 0; j < neqn; j++){
    qdots = DotProduct(de, &tgrad[j*3]);
    dq = (tvarsR[j] - tvarsL[j] - qdots)/ds2;
    for(i = 0; i < 3; i++){
      tgrad[j*3 + i] += dq*de[i];
    }
  }
    
  //call the diffusive member function for our current turbulence model
  turb->Diffusive(nu, tgrad, tvarsL, tvarsR, avec, dgrad, tresL, tresR, tjacL, tjacR);

  Type* resL = &turb->crs.b[left_cv*neqn];
  Type* resR = &turb->crs.b[right_cv*neqn];
  Type* jacL = turb->crs.A->GetPointer(right_cv, left_cv);
  Type* jacR = turb->crs.A->GetPointer(left_cv, right_cv);

  for(i = 0; i < neqn; i++){
    //positive on the left b/c of outward pointing normal and div. thm.
    resL[i] += tresL[i];
    resR[i] -= tresR[i];
  }

  for(i = 0; i < neqn*neqn; i++){
    jacL[i] -= tjacL[i];
    jacR[i] -= tjacR[i];
  }

  *size = 0;
  *ptrL = NULL;
  *ptrR = NULL;

  return;
}


template <class Type>
void Bkernel_Diffusive(B_KERNEL_ARGS)
{  
  Int i, j;
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  Int neqn = turb->neqn;

  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;

  Type* tvarsL = &turb->tvar[left_cv*neqn];
  Type* tvarsR = &turb->tvar[right_cv*neqn];

  Type* tresL = (Type*)alloca(sizeof(Type)*neqn);
  Type* tresR = (Type*)alloca(sizeof(Type)*neqn);
  Type* tjacL = (Type*)alloca(sizeof(Type)*neqn*neqn);
  Type* tjacR = (Type*)alloca(sizeof(Type)*neqn*neqn);

  Type dgrad;
  Type ds2;
  Type* de = (Type*)alloca(sizeof(Type)*3);
  Type* xL = &m->cg[left_cv*3];
  Type* xR = &m->cg[right_cv*3];

  ds2 = 0.0;
  for(i = 0; i < 3; i++){
    de[i] = xR[i] - xL[i];
    ds2 += de[i]*de[i];
  }
  
  Type dx, dy, dz, d;
  dx = de[0]*avec[0];
  dy = de[1]*avec[1];
  dz = de[2]*avec[2];

  d = dx + dy + dz;
  dgrad = d/ds2;

  Int qneqn = eqnset->neqn;
  Int qnvars = qneqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*qnvars];
  Type* qR = &space->q[right_cv*qnvars];
  Type* qavg = (Type*)alloca(sizeof(Type)*qnvars);

  for(i = 0; i < qneqn; i++){
    qavg[i] = 0.5*(qL[i] + qR[i]);
  }
  eqnset->ComputeAuxiliaryVariables(qavg);
  Type rho = eqnset->GetDensity(qavg);
  Type mu = eqnset->ComputeViscosity(qavg);
  Type nu = mu/rho;

  Type* tgrad = (Type*)alloca(sizeof(Type)*neqn*3);
  if(m->IsGhostNode(right_cv)){
    for(i = 0; i < neqn*3; i++){
      tgrad[i] = 0.5*(turb->tgrad[left_cv*neqn*3 + i] + turb->tgrad[right_cv*neqn*3 + i]);
    }
    //do the h.o. extrapolation
    Type qdots, dq;
    for(j = 0; j < neqn; j++){
      qdots = DotProduct(de, &tgrad[j*3]);
      dq = (tvarsR[j] - tvarsL[j] - qdots)/ds2;
      for(i = 0; i < 3; i++){
	tgrad[j*3 + i] += dq*de[i];
      }
    }
  }
  else{
    //if we are on a physical boundary, just use the grad computed at the wall
    //doing extrapolations has shown to be unstable
    for(i = 0; i < neqn*3; i++){
      tgrad[i] = turb->tgrad[left_cv*neqn*3 + i];
    }
  }

  //call the diffusive member function for our current turbulence model
  turb->Diffusive(nu, tgrad, tvarsL, tvarsR, avec, dgrad, tresL, tresR, tjacL, tjacR);

  //now put the returned residuals and jacobians in the correct place
  for(i = 0; i < neqn; i++){
    //positive on the left b/c of outward pointing normal and div. thm.
    turb->crs.b[left_cv*neqn + i] += tresL[i];
  }

  if(m->IsGhostNode(right_cv)){
    Type* jacR = turb->crs.A->GetPointer(left_cv, right_cv);
    for(i = 0; i < neqn*neqn; i++){
      jacR[i] -= tjacR[i];
    }
  }
  else{
    Type* jacL = turb->crs.A->GetPointer(left_cv, left_cv);
    for(i = 0; i < neqn*neqn; i++){
      jacL[i] += tjacL[i];
    }
  }

  *size = 0;
  *ptrL = NULL;

  return;
}

//pass through for type resolution
template <class Type>
void BC_Turb_Kernel(B_KERNEL_ARGS)
{
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  turb->BC_Kernel(space, cvid, left_cv, right_cv, avec, vdotn, velw, ptrL, ptrR, 
		  tempL, tempR, size, custom, eid, factag);
  return;
}

//pass through for type resolution
template <class Type>
void BC_Turb_Jac_Kernel(B_KERNEL_ARGS)
{
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  turb->BC_Jac_Kernel(space, cvid, left_cv, right_cv, avec, vdotn, velw, ptrL, ptrR, 
		      tempL, tempR, size, custom, eid, factag);
  return;
}
