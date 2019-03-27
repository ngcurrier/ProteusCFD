#include "eqnset.h"
#include "driver.h"
#include "solutionField.h"
#include "solutionSpace.h"
#include "bc.h"
#include "bcobj.h"
#include "param.h"
#include "mesh.h"
#include "matrix.h"
#include "residual.h"

template <class Type>
void ComputeJacobians(SolutionSpace<Type>* space)
{
  ComputeSpatialJacobian(space);
  //only add temporal terms if solution is unsteady
  ContributeTemporalTerms(space);

  return;
}

template <class Type>
void Compute_dRdQ(CRS<Type>& crs, SolutionSpace<Type>* space, Int type)
{
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* param = space->param;
  Int i, j, k;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nnode = m->GetNumNodes();

  //zero old jacobians for accumulation
  crs.BlankMatrix();

  Int temp = param->boundaryJacEval;
  //we have to use exact boundary jacobians to get dR/dQ
  param->boundaryJacEval = 0;

  void (*numjac) (KERNEL_ARGS);
  void (*boundary) (B_KERNEL_ARGS);

  if(type == 0){
    //upwind
    numjac = Kernel_NumJac;
    boundary = Bkernel_NumJac;
  }
  else if(type == 1){ 
    //central
    numjac = Kernel_NumJac_Centered;
    boundary = Bkernel_NumJac_Centered;
  }
  else if(type == 2){
    //complex
    numjac = Kernel_NumJac_Complex;
    boundary = Bkernel_NumJac_Complex;
  }
  else if(type == 3){
    numjac = Kernel_NumJac_Centered;
    boundary = Bkernel_NumJac_Complex;
  }
  else{
    std::cerr << "WARNING: invalid type selection for dR/dQ calculation" << std::endl;
    numjac = Kernel_NumJac_Centered;
    boundary = Bkernel_NumJac_Centered;
  }

  Kernel<Type> NumJac(numjac);
  Kernel<Type> Diag_NumJac(Kernel_Diag_NumJac);
  Kernel<Type> BNumJac(boundary);

  //compute off diagonal jacobians
  Driver(space, NumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);
  //contribute boundary jacobian parts to diagonal and compute the
  //off-diagonal jacobians for the parallel edges
  Bdriver(space, BNumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);
  //add in viscous terms if necessary
  if(param->viscous){
    Kernel<Type> Visc(Kernel_Viscous_Jac);
    Kernel<Type> BVisc(Bkernel_Viscous_Jac);
    Driver(space, Visc, eqnset->neqn*eqnset->neqn, (void*)&crs);
    Bdriver(space, BVisc, eqnset->neqn*eqnset->neqn, (void*)&crs);
  }
  //sum off diagonal jacobians to diagonal
  Driver(space, Diag_NumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);

  //compute source term jacobians
  RCmplx* qc = new RCmplx[nvars];
  RCmplx* source = new RCmplx[nvars];
  RCmplx h(0.0, 1.0e-11);
  EqnSet<RCmplx >* ceqnset = space->ceqnset;
  for(i = 0; i < nnode; i++){
    Type* jac = crs.A->GetPointer(i, i);
    for(j = 0; j < neqn; j++){
      for(k = 0; k < neqn; k++){
	qc[k] = space->q[i*nvars + k];
      }
      qc[j] += h;
      ceqnset->ComputeAuxiliaryVariables(qc);
      ceqnset->SourceTerm(qc, m->vol[i], source); 
      for(k = 0; k < neqn; k++){
	//this is negative for dR/dQ b/c the source term function 
	//itself is written as if on the RHS of the equations already
	jac[k*neqn + j] -= imag(source[k])/imag(h);
      }
    }
  }

  delete [] qc;
  delete [] source;

  //this is a hook which allow us to modify the boundary nodes' jacobians
  Kernel<Type> JacModifyBC(Bkernel_BC_Jac_Modify);
	  BdriverNoScatter(space, JacModifyBC, eqnset->neqn*eqnset->neqn, (void*)&crs);

  //set the flag back to what is necessary for solution
  eqnset->param->boundaryJacEval = temp;


  return;
}

template <class Type>
void Compute_dRdQ_Transpose(CRS<Type>& crs, SolutionSpace<Type>* space, Int type)
{
  //compute dRdQ according to type
  Compute_dRdQ(crs, space, type);

  //and compute the transpose
  crs.A->CRSTranspose();

  return;
}

template <class Type>
void ComputeSpatialJacobian(SolutionSpace<Type>* space)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* param = eqnset->param;
  CRS<Type>& crs = *space->crs;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nnode = m->GetNumNodes();

  //zero old jacobians for accumulation
  crs.BlankMatrix();

  void (*numjac) (KERNEL_ARGS);
  void (*bnumjac) (B_KERNEL_ARGS);

  if(param->fieldJacType == 0){
    numjac = Kernel_NumJac;
  }
  else if(param->fieldJacType == 1){
    numjac = Kernel_NumJac_Centered;
  }
  else if(param->fieldJacType == 2){
    numjac = Kernel_NumJac_Complex;
  }
  else{
    std::cerr << "WARNING: invalid type selection for field jacobian calculation" << std::endl;
    numjac = Kernel_NumJac;
  }
  
  //
  //PICK YOUR POISON FOR THE BOUNDARY CONTRIBUTIONS
  //
  if(param->boundaryJacType == 0){
    bnumjac = Bkernel_NumJac;
  }
  else if(param->boundaryJacType == 1){
    bnumjac = Bkernel_NumJac_Centered;
  }
  else if(param->boundaryJacType == 2){
    bnumjac = Bkernel_NumJac_Complex;
  }
  else{
    std::cerr << "WARNING: invalid type selection for field jacobian calculation" << std::endl;
    bnumjac = Bkernel_NumJac;
  }

  Kernel<Type> NumJac(numjac);
  Kernel<Type> Diag_NumJac(Kernel_Diag_NumJac);
  Kernel<Type> BNumJac(bnumjac);

  //compute off diagonal jacobians
  Driver(space, NumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);
  //contribute boundary jacobian parts to diagonal and compute the
  //off-diagonal jacobians for the parallel edges
  Bdriver(space, BNumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);
  //add in viscous terms if necessary
  if(param->viscous){
    Kernel<Type> Visc(Kernel_Viscous_Jac);
    Kernel<Type> BVisc(Bkernel_Viscous_Jac);
    Driver(space, Visc, eqnset->neqn*eqnset->neqn, (void*)&crs);
    Bdriver(space, BVisc, eqnset->neqn*eqnset->neqn, (void*)&crs);
  }
  //sum off diagonal jacobians to diagonal
  Driver(space, Diag_NumJac, eqnset->neqn*eqnset->neqn, (void*)&crs);

  //compute source term jacobians
  Type* A = new Type[neqn2];
  for(i = 0; i < nnode; i++){
    eqnset->SourceTermJacobian(&space->q[i*nvars], m->vol[i], A);
    Type* jac = crs.A->GetPointer(i, i);
    for(j = 0; j < neqn2; j++){
      //this is negative for dR/dQ b/c the source term function 
      //itself is written as if on the RHS of the equations already
      jac[j] -= A[j];
    }
  }
  delete [] A;
}


template <class Type>
void ContributeTemporalTerms(SolutionSpace<Type>* space)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* param = space->param; //add temporal diagonal contributions
  CRS<Type>& crs = *space->crs;
  Type* dt = space->GetFieldData("timestep", FIELDS::STATE_NONE);
  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type cnp1 = 1.0;
  Int torder = param->torder;
  Int nnode = m->GetNumNodes();

  if(space->iter > 1 &&  torder == 2){
    cnp1 = 1.5;
  }
#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(m, cnp1, neqn, eqnset, dt, param, space, nvars) \
  private(i)
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif    
    for(i = 0; i < nnode; i++){
      eqnset->ContributeTemporalTerms(&space->q[i*nvars], m->vol[i], cnp1, (Type)param->dt,
				      dt[i], crs.A->GetPointer(i,i), beta[i]);
    }
  }

  //this is a hook which allow us to modify the boundary nodes' jacobians
  Kernel<Type> JacModifyBC(Bkernel_BC_Jac_Modify);
  BdriverNoScatter(space, JacModifyBC, eqnset->neqn*eqnset->neqn, (void*)&crs);
}


template <class Type>
void Kernel_NumJac(KERNEL_ARGS){
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type* QPR = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPL = (Type*)alloca(sizeof(Type)*nvars);
  Type* fluxS = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxL = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxR = (Type*)alloca(sizeof(Type)*neqn);
  Type h = 1.0e-8;

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];
  Type betaR = beta[right_cv];
  Type avbeta = 0.5*(betaL + betaR);

  //get reference state
  eqnset->NumericalFlux(QL, QR, avec, vdotn, fluxS, avbeta);
  
  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);

    //perturb variables
    QPL[i] += h;
    QPR[i] += h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->NumericalFlux(QPL, QR, avec, vdotn, fluxL, avbeta);
    eqnset->NumericalFlux(QL, QPR, avec, vdotn, fluxR, avbeta);
          
    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = real(fluxS[j] - fluxL[j])/h;
    }
    for(j = 0; j < neqn; j++){
      tempR[j*neqn + i] = real(fluxR[j] - fluxS[j])/h;
    }
  }

  //set data necessary for driver scatter
  *size = neqn2;

  *ptrR = crs.A->GetPointer(left_cv, right_cv);
  *ptrL = crs.A->GetPointer(right_cv, left_cv);

  return;
}


template <class Type>
void Kernel_NumJac_Centered(KERNEL_ARGS){
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type* QPR = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPL = (Type*)alloca(sizeof(Type)*nvars);
  Type* fluxLd = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxRd = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxLu = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxRu = (Type*)alloca(sizeof(Type)*neqn);
  Type h = 1.0e-8;

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];
  Type betaR = beta[right_cv];
  Type avbeta = 0.5*(betaL + betaR);

  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);

    //perturb variables
    QPL[i] += h;
    QPR[i] += h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->NumericalFlux(QPL, QR, avec, vdotn, fluxLu, avbeta);
    eqnset->NumericalFlux(QL, QPR, avec, vdotn, fluxRu, avbeta);

    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);

    //perturb variables
    QPL[i] -= h;
    QPR[i] -= h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->NumericalFlux(QPL, QR, avec, vdotn, fluxLd, avbeta);
    eqnset->NumericalFlux(QL, QPR, avec, vdotn, fluxRd, avbeta);
          
    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = (fluxLd[j] - fluxLu[j])/(2.0*h);
    }
    for(j = 0; j < neqn; j++){
      tempR[j*neqn + i] = (fluxRu[j] - fluxRd[j])/(2.0*h);
    }
  }

  //set data necessary for driver scatter
  *size = neqn2;
  *ptrR = crs.A->GetPointer(left_cv, right_cv);
  *ptrL = crs.A->GetPointer(right_cv, left_cv);

  return;
}

template <class Type>
void Kernel_NumJac_Complex(KERNEL_ARGS){
  Int i, j;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  CRS<Type>& crs = *(CRS<Type>*)custom;
  RCmplx* QR = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QL = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QPR = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QPL = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* fluxL = (RCmplx*)alloca(sizeof(RCmplx)*neqn);
  RCmplx* fluxR = (RCmplx*)alloca(sizeof(RCmplx)*neqn);
  RCmplx h (0.0, 1.0e-11);
  RCmplx avecC[4];
  RCmplx vdotnC = vdotn;

  //get complex eqnset from param creation
  EqnSet<RCmplx>* ceqnset = space->ceqnset;

  for(j = 0; j < nvars; j++){
    QR[j] = qR[j];
    QL[j] = qL[j];
  }
  for(j = 0; j < 4; j++){
    avecC[j] = avec[j];
  }

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];
  Type betaR = beta[right_cv];
  Type avbeta = 0.5*(betaL + betaR);

  RCmplx cavbeta = avbeta;
    
  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(RCmplx)*nvars);
    memcpy(QPR, QR, sizeof(RCmplx)*nvars);

    //perturb variables
    QPL[i] += h;
    QPR[i] += h;
    ceqnset->ComputeAuxiliaryVariables(QPL);
    ceqnset->ComputeAuxiliaryVariables(QPR);

    ceqnset->NumericalFlux(QPL, QR, avecC, vdotnC, fluxL, cavbeta);
    ceqnset->NumericalFlux(QL, QPR, avecC, vdotnC, fluxR, cavbeta);

    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = -imag(fluxL[j])/imag(h);
    }
    for(j = 0; j < neqn; j++){
      tempR[j*neqn + i] = imag(fluxR[j])/imag(h);
    }
  }

  //set data necessary for driver scatter
  *size = neqn2;
  *ptrR = crs.A->GetPointer(left_cv, right_cv);
  *ptrL = crs.A->GetPointer(right_cv, left_cv);

  return;
}


template <class Type>
void Kernel_Diag_NumJac(KERNEL_ARGS){
  Int i;
  CRS<Type>* crs = (CRS<Type>*) custom;
  Int neqn = crs->neqn;
  Int neqn2 = neqn*neqn;

  //set data necessary for driver scatter
  *size = neqn2;
  *ptrL = crs->A->GetPointer(left_cv, left_cv);
  *ptrR = crs->A->GetPointer(right_cv, right_cv);

  Type* jacL = crs->A->GetPointer(right_cv, left_cv);
  Type* jacR = crs->A->GetPointer(left_cv, right_cv);

  //sum the off-diagonals on the column up to the diagonal
  for(i = 0; i < neqn2; i++){
    tempL[i] = -jacL[i];
    tempR[i] = -jacR[i];
  }
  return;
}

template <class Type>
void Bkernel_NumJac(B_KERNEL_ARGS){
  Int i, j;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type* QPR = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPL = (Type*)alloca(sizeof(Type)*nvars);
  Type* fluxS = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxL = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxR = (Type*)alloca(sizeof(Type)*neqn);
  Type h = 1.0e-8;

  BoundaryConditions<Real>* bc = space->bc;
  Int bcType = bc->GetBCType(factag);
  BCObj<Real>* bcobj = bc->GetBCObj(factag);
  Int bcId; 

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];

  Type* Qref = (Type*)alloca(sizeof(Type)*nvars);
  bcobj->GetQref(Qref);
  CalculateBoundaryVariables(eqnset, m, space, QL, QR, Qref, avec, bcType, eid, bcobj, vdotn, velw, factag);
  //get reference state
  eqnset->BoundaryFlux(QL, QR, avec, vdotn, fluxS, betaL);
  
  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);
    
    //perturb internal variables
    QPL[i] += h;
    QPR[i] += h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->BoundaryFlux(QL, QPR, avec, vdotn, fluxR, betaL);

    if(!param->boundaryJacEval && !m->IsGhostNode(right_cv)){
      memcpy(QPR, QR, sizeof(Type)*nvars);
      eqnset->ComputeAuxiliaryVariables(QPR);
      CalculateBoundaryVariables(eqnset, m, space, QPL, QPR, Qref, avec, bcType, eid, bcobj, vdotn, velw, factag);
      eqnset->BoundaryFlux(QPL, QPR, avec, vdotn, fluxL, betaL);
    }
    else{
      eqnset->BoundaryFlux(QPL, QR, avec, vdotn, fluxL, betaL);    
    }
           
    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = (fluxL[j] - fluxS[j])/h;
      tempR[j*neqn + i] = (fluxR[j] - fluxS[j])/h;
    }
  }
  
  //sum into boundary jacobian spot
  if(m->IsGhostNode(right_cv)){
    //this is the standard contribution
    *size = neqn2;
    *ptrR = crs.A->GetPointer(cvid, right_cv);

    //this is a contribution to the diagonal to avoid a parallel
    //communication step
    *ptrL = crs.A->GetPointer(left_cv, left_cv);
    return;
  }

  //set data necessary for driver scatter to diagonal jacobian
  *size = neqn2;
  *ptrL = crs.A->GetPointer(cvid, cvid);
  *ptrR = NULL;
  return;
}

template <class Type>
void Bkernel_NumJac_Centered(B_KERNEL_ARGS){
  Int i, j;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type* QPR = (Type*)alloca(sizeof(Type)*nvars);
  Type* QPL = (Type*)alloca(sizeof(Type)*nvars);
  Type* fluxLu = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxRu = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxLd = (Type*)alloca(sizeof(Type)*neqn);
  Type* fluxRd = (Type*)alloca(sizeof(Type)*neqn);
  Type h = 1.0e-8;

  BoundaryConditions<Real>* bc = space->bc;
  BCObj<Real>* bcobj = bc->GetBCObj(factag);
  Int bcType = bc->GetBCType(factag); 
  Type* Qref = (Type*)alloca(sizeof(Type)*nvars);
  bcobj->GetQref(Qref);

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];

  CalculateBoundaryVariables(eqnset, m, space, QL, QR, Qref, avec, bcType, eid, bcobj, vdotn, velw, factag);

  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);
    
    //perturb internal variables
    QPL[i] += h;
    QPR[i] += h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->BoundaryFlux(QL, QPR, avec, vdotn, fluxRu, betaL);

    if(!param->boundaryJacEval && !m->IsGhostNode(right_cv)){
      memcpy(QPR, QR, sizeof(Type)*nvars);
      eqnset->ComputeAuxiliaryVariables(QPR);
      //calculate new external variables using
      //applied BCs
      CalculateBoundaryVariables(eqnset, m, space, QPL, QPR, Qref, avec, bcType, eid, bcobj, vdotn, velw, factag);
      eqnset->BoundaryFlux(QPL, QPR, avec, vdotn, fluxLu, betaL);
    }
    else{
      eqnset->BoundaryFlux(QPL, QR, avec, vdotn, fluxLu, betaL);    
    }

    memcpy(QPL, QL, sizeof(Type)*nvars);
    memcpy(QPR, QR, sizeof(Type)*nvars);
    
    //perturb internal variables
    QPL[i] -= h;
    QPR[i] -= h;
    eqnset->ComputeAuxiliaryVariables(QPL);
    eqnset->ComputeAuxiliaryVariables(QPR);

    eqnset->BoundaryFlux(QL, QPR, avec, vdotn, fluxRd, betaL);
    
    if(param->boundaryJacEval && !m->IsGhostNode(right_cv)){
      memcpy(QPR, QR, sizeof(Type)*nvars);
      eqnset->ComputeAuxiliaryVariables(QPR);
      //calculate new external variables using
      //applied BCs
      CalculateBoundaryVariables(eqnset, m, space, QPL, QPR, Qref, avec, bcType, eid, bcobj, vdotn, velw, factag);
      eqnset->BoundaryFlux(QPL, QPR, avec, vdotn, fluxLd, betaL);
    }
    else{
      eqnset->BoundaryFlux(QPL, QR, avec, vdotn, fluxLd, betaL);    
    }    
          
    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = (fluxLu[j] - fluxLd[j])/(2.0*h);
      tempR[j*neqn + i] = (fluxRu[j] - fluxRd[j])/(2.0*h);
    }
  }

  //sum into boundary jacobian spot
  if(m->IsGhostNode(right_cv)){
    //this is the standard contribution
    *size = neqn2;
    *ptrR = crs.A->GetPointer(cvid, right_cv);

    //this is a contribution to the diagonal to avoid a parallel
    //communication step
    *ptrL = crs.A->GetPointer(left_cv, left_cv);
    return;
  }
  
  //set data necessary for driver scatter to diagonal jacobian
  *size = neqn2;
  *ptrL = crs.A->GetPointer(cvid, cvid);
  *ptrR = NULL;
 
  return;
}

template <class Type>
void Bkernel_NumJac_Complex(B_KERNEL_ARGS){
  Int i, j;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  RCmplx* QR = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QL = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QPR = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* QPL = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  RCmplx* fluxL = (RCmplx*)alloca(sizeof(RCmplx)*neqn);
  RCmplx* fluxR = (RCmplx*)alloca(sizeof(RCmplx)*neqn);
  RCmplx avecC[4];
  RCmplx vdotnC = vdotn;
  RCmplx h(0.0, 1.0e-11);

  //get complex eqnset from param creation
  EqnSet<RCmplx>* ceqnset = space->ceqnset;

  BoundaryConditions<Real>* bc = space->bc;
  BCObj<Real>* bcobj = bc->GetBCObj(factag);
  Int bcType = bc->GetBCType(factag); 

  RCmplx* Qref = (RCmplx*)alloca(sizeof(RCmplx)*nvars);
  bcobj->GetQref(Qref);

  Type* beta = space->GetFieldData("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];
  RCmplx cbetaL = betaL;

  for(j = 0; j < 4; j++){
    avecC[j] = avec[j];
  }
  for(j = 0; j < nvars; j++){
    QL[j] = qL[j];
    QR[j] = qR[j];
  }

  CalculateBoundaryVariables(ceqnset, m, space, QL, QR, Qref, avecC, bcType, eid, bcobj, vdotn, velw, factag);
  
  for(i = 0; i < neqn; i++){
    memcpy(QPL, QL, sizeof(RCmplx)*nvars);
    memcpy(QPR, QR, sizeof(RCmplx)*nvars);
    
    //perturb internal variables
    QPL[i] += h;
    QPR[i] += h;
    ceqnset->ComputeAuxiliaryVariables(QPL);
    ceqnset->ComputeAuxiliaryVariables(QPR);

    ceqnset->BoundaryFlux(QL, QPR, avecC, vdotnC, fluxR, cbetaL);

    if(!param->boundaryJacEval && !m->IsGhostNode(right_cv)){
      memcpy(QPR, QR, sizeof(RCmplx)*nvars);
      ceqnset->ComputeAuxiliaryVariables(QPR);
      //calculate new external variables using applied BCs
      CalculateBoundaryVariables(ceqnset, m, space, QPL, QPR, Qref, avecC, bcType, eid, bcobj, vdotn, velw, factag);
      ceqnset->BoundaryFlux(QPL, QPR, avecC, vdotnC, fluxL, cbetaL);
    }
    else{
      ceqnset->BoundaryFlux(QPL, QR, avecC, vdotnC, fluxL, cbetaL);    
    }
    
    for(j = 0; j < neqn; j++){
      tempL[j*neqn + i] = imag(fluxL[j])/imag(h);
      tempR[j*neqn + i] = imag(fluxR[j])/imag(h);
    }
  }

  //sum into boundary jacobian spot
  if(m->IsGhostNode(right_cv)){
    //this is the standard contribution
    *size = neqn2;
    *ptrR = crs.A->GetPointer(cvid, right_cv);

    //this is a contribution to the diagonal to avoid a parallel
    //communication step
    *ptrL = crs.A->GetPointer(left_cv, left_cv);
    
    return;
  }
  
  //set data necessary for driver scatter to diagonal jacobian
  *size = neqn2;
  *ptrL = crs.A->GetPointer(cvid, cvid);
  *ptrR = NULL;
}


template <class Type>
void Kernel_Viscous_Jac(KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Type* mut = space->GetFieldData("mut", FIELDS::STATE_NONE);
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* xL;
  Type* xR;
  Type* dx = (Type*)alloca(sizeof(Type)*3);

  xL = m->xyz + 3*left_cv;
  xR = m->xyz + 3*right_cv;

  //get turbulent viscosity
  Type tmut = (mut[left_cv] + mut[right_cv])/2.0;
  
  Type s2 = 0.0;
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }

  //call viscous flux routine
  eqnset->ViscousJacobian(qL, qR, dx, s2, avec, tmut, tempL, tempR);

  //set data necessary for driver scatter
  *size = neqn2;
  *ptrL = crs.A->GetPointer(right_cv, left_cv);
  *ptrR = crs.A->GetPointer(left_cv, right_cv);
}


template <class Type>
void Bkernel_Viscous_Jac(B_KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  CRS<Type>& crs = *(CRS<Type>*)custom;
  Type* mut = space->GetFieldData("mut", FIELDS::STATE_NONE);
  Int neqn = eqnset->neqn;
  Int neqn2 = neqn*neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* xL;
  Type* xR;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type s2 = 0.0;

  BoundaryConditions<Real>* bc = space->bc;
  Int bcType = bc->GetBCType(factag);

  xL = m->xyz + 3*left_cv;
  xR = m->xyz + 3*right_cv;
  
  //get turbulent viscosity
  Type tmut;
  if(m->IsGhostNode(right_cv)){
    tmut = (mut[left_cv] + mut[right_cv])/2.0;
  }
  else{
   tmut = mut[left_cv];
  }

  //if we have gradient info on a parallel interface
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }

  //call viscous flux routine
  eqnset->ViscousJacobian(qL, qR, dx, s2, avec, tmut, tempL, tempR);

  //sum into boundary jacobian spot
  if(m->IsGhostNode(right_cv)){
    //this is the standard contribution
    *size = neqn2;
    *ptrR = crs.A->GetPointer(cvid, right_cv);

    //this is a contribution to the diagonal to avoid a parallel
    //communication step
    for(i = 0; i < neqn2; i++){
      tempL[i] = -tempL[i];
    }
    *ptrL = crs.A->GetPointer(left_cv, left_cv);
    
    return;
  }

  //set data necessary for driver scatter
  *size = neqn2;
  *ptrL = crs.A->GetPointer(left_cv, left_cv);

  //correct signs (tempspace)
  for(i = 0; i < neqn*neqn; i++){
    tempL[i] = -tempL[i];
  }

  return;
}

#if 0
//This function will compute the diagonal jacobian for point (ptid) using
//the full blown residual call SpatialResidual() and CTSE, then places it in JacCheck
//This is super expensive and should only be used to sanity check dR/dQ routines
template <class Type>
void Diag_Jacobian_Check(Real* JacCheck, SolutionSpace<Type>* space, Int ptid){

  Int i,j;
  RCmplx perturb;
  perturb.real() = 0.0;
  perturb.imag() = 1.0e-8;

  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  
  SolutionSpace<RCmplx> cSpace(*space);
  //NOTE: no need to set parallel pointers here since residual call is local only

  for(j = 0; j < neqn; j++){
    cSpace.q[ptid*nvars + j] += perturb;
    //compute spatial residual
    UpdateBCs(&cSpace);
    SpatialResidual(&cSpace);
    for(i = 0; i < neqn; i++){
      JacCheck[i*neqn + j] = imag(cSpace.crs->b[ptid*neqn + i])/imag(perturb);
    }
    cSpace.q[ptid*nvars + j] -= perturb;
  }
 

  delete ceqnset;

  return;
}
#endif
