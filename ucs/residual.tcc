#include "bc.h"
#include "limiters.h"
#include "solutionSpace.h"
#include "solutionField.h"
#include "driver.h"
#include "eqnset.h"
#include "param.h"
#include "mesh.h"
#include "gradient.h"
#include "gaussian.h"

template <class Type>
std::vector<Type> ComputeResiduals(SolutionSpace<Type>* space)
{
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;

  //compute spatial residual
  SpatialResidual(space);

  //if solution method is implicit add in temporal residual
  //also only add temporal terms if solution is unsteady
  if(param->torder){
    TemporalResidual(space, (Type)param->dt, space->iter, param->torder);
    if(param->gcl){
      if(param->gcl == 1){
	//this version uses the integrated face speeds
	ContributeGCL(space, (Type)param->dt, space->iter, param->torder);
      }
      else{
	//this version follows the strict derivation for a BDF2 dv/dt term
	ContributeGCL2(space, (Type)param->dt, space->iter, param->torder);
      }
    }
  }

  //this is a hook to modify any residuals directly due to boundary
  //conditions, most notably hardset viscous BCs
  Kernel<Type> ResModifyBC(Bkernel_BC_Res_Modify);
  BdriverNoScatter(space, ResModifyBC, eqnset->neqn, NULL);

  //compute residual norm before static sources, this gives a better indication
  //of whether or not the system is converged
  Type residGlobal = ParallelL2Norm(space->p, space->crs->b, space->m->nnode*neqn);

  std::vector<Type> res = StridedParallelL2Norm(space->p, space->crs->b, space->m->nnode, neqn);

  //this applies any constant source terms which do not vanish as the solution
  //converges. Just keeps their contributions out of the residual computation
  ExtraSourceResidual(space);

  std::vector<Type> resComponents(res.size() + 1);

  resComponents[0] = residGlobal;
  for(Int i = 0; i < neqn; i++){
    resComponents[1+i] = res[i];
  }

  return resComponents;
}

template <class Type>
void SpatialResidual(SolutionSpace<Type>* space)
{
  Int i, j;
  Mesh<Type>* m = space->m;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* source = new Type[neqn];
  Type sourcenorm = 0.0;

  //zero residuals
  space->crs->BlankB();
  
  Kernel<Type> Inviscid_Flux(Kernel_Inviscid_Flux);
  Kernel<Type> BInviscid_Flux(Bkernel_Inviscid_Flux);

  //call drivers to calculate inviscid residuals
  Driver(space, Inviscid_Flux, neqn*neqn, NULL);
  Bdriver(space, BInviscid_Flux, neqn*neqn, NULL);

  if(param->viscous){
    //pre-fetch "mut" to avoid the lookup at every node
    Type* mut = space->GetField("mut", FIELDS::STATE_NONE);

    Kernel<Type> Viscous_Flux(Kernel_Viscous_Flux);
    Kernel<Type> BViscous_Flux(Bkernel_Viscous_Flux);

    //call drivers to calculate viscous residuals
    Driver(space, Viscous_Flux, neqn*neqn, (void*)mut);
    Bdriver(space, BViscous_Flux, neqn*neqn, (void*)mut);

    if(param->errorTransport){
      space->GetField("ETEVisc").Fill(0.0);
      Kernel<Type> Viscous_Flux_ETE_Src(Kernel_Viscous_Src);
      Kernel<Type> BViscous_Flux_ETE_Src(Bkernel_Viscous_Src);
      Driver(space, Viscous_Flux_ETE_Src, neqn*neqn, (void*)mut);
      Bdriver(space, BViscous_Flux_ETE_Src, neqn*neqn, (void*)mut);
    }
  }

  //contribute source terms
  for(i = 0; i < m->nnode; i++){
    eqnset->SourceTerm(&space->q[i*nvars], m->vol[i], source);
    for(j = 0; j < neqn; j++){
      space->crs->b[i*neqn + j] += source[j];
      sourcenorm += source[j]*source[j];
    }
  }
  sourcenorm = sqrt(sourcenorm);
  if(sourcenorm != 0.0){
    std::cout << " ||Source||: " << sourcenorm << " ";
  }

  delete [] source;

  return;
}

template <class Type>
void TemporalResidual(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder)
{
  Int i, j, offset;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type cnp1, cnm1;
  Type dt, dtm1;
  Type* q = new Type[nvars];
  Type* dq = new Type[neqn];
  Type* dqm1 = new Type[neqn];

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
  for(i = 0; i < m->nnode; i++){
    dt = cnp1*m->vol[i]/timestep;           //use v_n+1
    if(eqnset->param->gcl && iter > 2){
      dtm1 = cnm1*m->vololdm1[i]/timestep;  //use v_n-1
    }
    else{
      dtm1 = cnm1*m->vol[i]/timestep;       //use v_n+1
    }
    offset = i*nvars;
    
    memcpy(q, &space->q[offset + 0], sizeof(Type)*nvars);
    //since we are contributing directly to the residual, we must make
    //sure that q is the actual q we write the equations for, i.e. conservative
    eqnset->NativeToConservative(q);
    for(j = 0; j < neqn; j++){
      //it is assumed here that the q's stored are conservative variables
      //and can be contributed directly to the RHS
      dq[j] = q[j] - space->qold[offset + j];
      dqm1[j] = space->qold[offset + j] - space->qoldm1[offset + j];
    }
    for(j = 0; j < neqn; j++){
      space->crs->b[i*neqn + j] -= dt*dq[j];
      space->crs->b[i*neqn + j] -= dtm1*dqm1[j];
    }
  }
    
  delete [] q;
  delete [] dq;
  delete [] dqm1;

  return;
}

template <class Type>
void ExtraSourceResidual(SolutionSpace<Type>* space)
{
  Param<Type>* param = space->param;

  //contribute any gaussian mass sources
  if(param->gaussianSource){
    space->gaussian->ApplyToResidual();
  }
}
template <class Type>
void Kernel_Inviscid_Flux(KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* param = eqnset->param;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* dQ = (Type*)alloca(sizeof(Type)*nvars);
  Type* xR,* xL;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* limiterL = &space->limiter->l[left_cv*neqn];
  Type* limiterR = &space->limiter->l[right_cv*neqn];
  Type chi = param->chi;

  Type* beta = space->GetField("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];
  Type betaR = beta[right_cv];
  Type avbeta = 0.5*(betaL + betaR);

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;
  
  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);

  if(param->sorder > 1 && (space->iter > param->nFirstOrderSteps)){
    Type* gradL = space->qgrad + 3*nterms*left_cv;
    Type* gradR = space->qgrad + 3*nterms*right_cv;

    //TODO: generalize incase we want to use a non-midpoint edge CV
    dx[0] = 0.5*(xR[0] - xL[0]);
    dx[1] = 0.5*(xR[1] - xL[1]);
    dx[2] = 0.5*(xR[2] - xL[2]);

    for(j = 0; j < neqn; j++){
      dQ[j] = qR[j] - qL[j];
    }

    eqnset->ExtrapolateVariables(QL, qL, dQ, gradL, dx, limiterL);

    for(j = 0; j < neqn; j++){
      dQ[j] = -dQ[j];
    }

    dx[0] = -dx[0];
    dx[1] = -dx[1];
    dx[2] = -dx[2];

    eqnset->ExtrapolateVariables(QR, qR, dQ, gradR, dx, limiterR);

    eqnset->ExtrapolatedToNative(QL);
    eqnset->ExtrapolatedToNative(QR);
    eqnset->ComputeAuxiliaryVariables(QL);
    eqnset->ComputeAuxiliaryVariables(QR);

#if 0
    //some debugging output for checking extrapolations
    std::cout << "CV_left: " << left_cv << " CV_right: " << right_cv 
	      << std::endl << std::endl;
    if(!m->IsInteriorNode(left_cv) || !m->IsInteriorNode(right_cv))
      {
	std::cout << "Other node found in interior kernel" << std::endl;
      };
    for(j = 0; j < neqn; j++){
      std::cout << "orig L: " << qL[j] << "  extrap L: " << QL[j] << std::endl;
      std::cout << "gradL " << gradL[j*3 + 0] << " " << gradL[j*3 + 1] << " " 
		<< gradL[j*3 + 2] << std::endl;
      std::cout << "orig R: " << qR[j] << "  extrap R: " << QR[j] << std::endl;
      std::cout << "gradR " << gradR[j*3 + 0] << " " << gradR[j*3 + 1] << " " 
		<< gradR[j*3 + 2] << std::endl << std::endl;
      std::cout << "limL: " << limiterL[j] << " limR: " << limiterR[j] << std::endl;
    }
#endif
    
  }

  eqnset->NumericalFlux(QL, QR, avec, vdotn, tempR, avbeta);

  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  *ptrR = &space->crs->b[right_cv*neqn];
  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempR[i];
  }

  return;
}

template <class Type>
void Bkernel_Inviscid_Flux(B_KERNEL_ARGS)
{
  Int i,j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Param<Type>* param = eqnset->param;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars); 
  Type* dQ = (Type*)alloca(sizeof(Type)*nvars);
 
  Type* xL;
  Type* xR;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type chi = param->chi;

  Type* beta = space->GetField("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];

  BoundaryConditions<Real>* bc = space->bc;
  Int bcNum = bc->bc_map[factag];
  Int bcId; 
  if(factag == 0){
    bcId = ParallelBoundary;
  }
  else{
    bcId = bc->bc_applied[bcNum];
  }  

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;
  
  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);

  if(param->sorder > 1 && (space->iter > param->nFirstOrderSteps)){
    
    //only do full right extrapolation if boundary node is ghost(parallel)
    if(m->IsGhostNode(right_cv)){
      Type* limiterL = &space->limiter->l[left_cv*neqn];
      Type* gradL = space->qgrad + 3*nterms*left_cv;
      
      for(j = 0; j < neqn; j++){
	dQ[j] = qR[j] - qL[j];
      }

      //TODO: generalize incase we want to use a non-midpoint edge CV
      dx[0] = 0.5*(xR[0] - xL[0]);
      dx[1] = 0.5*(xR[1] - xL[1]);
      dx[2] = 0.5*(xR[2] - xL[2]);

      eqnset->ExtrapolateVariables(QL, qL, dQ, gradL, dx, limiterL);
      
      for(j = 0; j < neqn; j++){
	dQ[j] = -dQ[j];
      }
      dx[0] = -dx[0];
      dx[1] = -dx[1];
      dx[2] = -dx[2];
      Type* gradR = space->qgrad + 3*nterms*right_cv;
      Type* limiterR = &space->limiter->l[right_cv*neqn];

      eqnset->ExtrapolateVariables(QR, qR, dQ, gradR, dx, limiterR);

    }
    eqnset->ExtrapolatedToNative(QL);
    eqnset->ExtrapolatedToNative(QR);
    eqnset->ComputeAuxiliaryVariables(QL);
    eqnset->ComputeAuxiliaryVariables(QR);
    
  }
  
  eqnset->BoundaryFlux(QL, QR, avec, vdotn, tempL, betaL);
  
  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempL[i];
  }

  return;
}

template <class Type>
void Kernel_Viscous_Flux(KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Type* mut = (Type*)custom;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* Qavg = (Type*)alloca(sizeof(Type)*nvars);
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* xR,* xL;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* gradL = space->qgrad + 3*nterms*left_cv;
  Type* gradR = space->qgrad + 3*nterms*right_cv;
  Type* grad = (Type*)alloca(sizeof(Type)*nterms*3);

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;

  //do dumb averaging for q vector... this affects the velocity terms internally
  //TODO: use h.o. extrapolations here
  for(i = 0; i < neqn; i++){
    Qavg[i] = (qL[i] + qR[i])/2.0;
  }
  eqnset->ComputeAuxiliaryVariables(Qavg);

  std::vector<Int> gradientsLoc;
  eqnset->GetGradientsLocation(gradientsLoc);
  if(gradientsLoc.size() != (UInt)nterms){
    Abort << "WARNING: gradients location vector does not match number of terms in gradient";
  }

  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);
  eqnset->NativeToExtrapolated(QL);
  eqnset->NativeToExtrapolated(QR);

  //get turbulent viscosity
  Type tmut = 0.5*(mut[left_cv] + mut[right_cv]);

  //use directional derivatives to get gradient at face
  //TODO: generalize in case we want to use a non-midpoint edge CV
  for(i = 0; i < nterms*3; i++){
    grad[i] = 0.5*(gradL[i] + gradR[i]);
  }

  Type s2 = 0.0;
  Type qdots, dq;
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }
  for(j = 0; j < nterms; j++){
    qdots = DotProduct(dx, &grad[j*3]);
    Int loc = gradientsLoc[j];
    dq = (QR[loc] - QL[loc] - qdots)/s2;
    for(i = 0; i < 3; i++){
      grad[j*3 + i] += dq*dx[i];
    }
  }

  //call viscous flux routine
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempR);


  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  *ptrR = &space->crs->b[right_cv*neqn];
  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempR[i];
  }

  return;
}

template <class Type>
void Bkernel_Viscous_Flux(B_KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Type* mut = (Type*)custom;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* Qavg = (Type*)alloca(sizeof(Type)*nvars);
  Type* xL;
  Type* xR;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* gradL = space->qgrad + 3*nterms*left_cv;
  Type* grad = (Type*)alloca(sizeof(Type)*nterms*3);

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;

  BoundaryConditions<Real>* bc = space->bc;
  Int bcNum = bc->bc_map[factag];
  Int bcId; 
  if(factag == 0){
    bcId = ParallelBoundary;
  }
  else{
    bcId = bc->bc_applied[bcNum];
  }  

  //use directional derivatives to get gradient at face
  //TODO: generalize in case we want to use a non-midpoint edge CV
  Type s2 = 0.0;
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }

  //do dumb averaging for q vector... this affects the velocity terms internally
  //TODO: use h.o. extrapolations here
  for(i = 0; i < neqn; i++){
    Qavg[i] = 0.5*(qL[i] + qR[i]);
  }
  eqnset->ComputeAuxiliaryVariables(Qavg);

  std::vector<Int> gradientsLoc;
  eqnset->GetGradientsLocation(gradientsLoc);
  if(gradientsLoc.size() != (UInt)nterms){
    Abort << "WARNING: gradients location vector does not match number of terms in gradient";
  }

  //get turbulent viscosity
  Type tmut;
  if(m->IsGhostNode(right_cv)){
    memcpy(QL, qL, sizeof(Type)*nvars);
    memcpy(QR, qR, sizeof(Type)*nvars);
    eqnset->NativeToExtrapolated(QL);
    eqnset->NativeToExtrapolated(QR);

    Type* gradR = space->qgrad + 3*nterms*right_cv;
    tmut = (mut[left_cv] + mut[right_cv])/2.0;
    for(i = 0; i < nterms*3; i++){
      grad[i] = 0.5*(gradL[i] + gradR[i]);
    }
    Type qdots, dq;
    for(j = 0; j < nterms; j++){
      qdots = DotProduct(dx, &grad[j*3]);
      Int loc = gradientsLoc[j];
      dq = (QR[loc] - QL[loc] - qdots)/s2;
      for(i = 0; i < 3; i++){
	grad[j*3 + i] += dq*dx[i];
      }
    }
  }
  else{
    //if we are on a boundary condition, this needs to be handled via a directional
    //derivative type rule, since our gradients do not include the boundary points
    tmut = mut[left_cv];
    //DirectionalDerivativeGrad(dx, s2, avec, qL, qR, grad, list, nterms);
    memcpy(grad, gradL, sizeof(Type)*3*nterms);
  }

  //call viscous flux routine
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempL);

  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempL[i];
  }

  return;
}

template <class Type>
void ContributeGCL2(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder)
{
  Int i, j;
  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int offset;
  Type cnp1, cnm1, cn;
  Type K;
  Type* q = new Type[neqn];

  Type res;
  Type maxres;

  if(iter > 2 && torder == 2){
    //coefficients for BDF2 are phi_n+1 = 1.5, phi_n = -2.0, phi_n-1 = 0.5 
    cnp1 = 1.5;
    cn = -2.0;
    cnm1 = -0.5;
  }
  else{
    cnp1 = 1.0;
    cn = -1.0;
    cnm1 = 0.0;
  }

  res = 0.0;
  maxres = 0.0;
  for(i = 0; i < m->nnode; i++){
    offset = i*nvars;
    K = (cnp1*m->vol[i] + cn*m->volold[i] - cnm1*m->vololdm1[i])/timestep;

    memcpy(q, &space->qold[offset + 0], sizeof(Type)*neqn);
    //since we are contributing directly to the residual, we must make
    //sure that q is the actual q we write the equations for, i.e. conservative
    for(j = 0; j < neqn; j++){
      space->crs->b[i*neqn + j] -= K*q[j];
      res += K*q[j]*K*q[j];
      maxres = MAX(K*q[j], maxres);
    }
  }
  
  res = sqrt(res)/(Type)m->nnode;
  std::cout << "\t||GCL-res||: " << res; 
  std::cout << "\tmax||GCL-res||: " << maxres;
  std::cout << "\t";

  delete [] q;
  
  return;
}

template <class Type>
void ContributeGCL(SolutionSpace<Type>* space, Type timestep, Int iter, Int torder)
{
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;

  Kernel<Type> GCL(GCL_Kernel);
  Kernel<Type> BGCL(GCL_Bkernel);
  Driver(space, GCL, neqn, NULL);
  Bdriver(space, BGCL, neqn, NULL);
  return;
}

template <class Type>
void GCL_Kernel(KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->qold[left_cv*nvars];
  Type* QR = &space->qold[right_cv*nvars];

  for(i = 0; i < neqn; i++){
    tempL[i] = avec[3]*vdotn*QL[i];
    tempR[i] = -avec[3]*vdotn*QR[i];
  }
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  *ptrR = &space->crs->b[right_cv*neqn];

  return;
}

template <class Type>
void GCL_Bkernel(B_KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->qold[left_cv*nvars];
  Type* QR = &space->qold[right_cv*nvars];
  
  for(i = 0; i < neqn; i++){
    tempL[i] = avec[3]*vdotn*QL[i];
    tempR[i] = -avec[3]*vdotn*QR[i];
  }
  *size = neqn;
  *ptrL = &space->crs->b[left_cv*neqn];
  *ptrR = NULL;

  return;
}

template <class Type>
void Kernel_Viscous_Src(KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Type* mut = (Type*)custom;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* Qavg = (Type*)alloca(sizeof(Type)*nvars);
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* xR,* xL;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* gradL = space->qgrad + 3*nterms*left_cv;
  Type* gradR = space->qgrad + 3*nterms*right_cv;
  Type* grad = (Type*)alloca(sizeof(Type)*nterms*3);

  Type* tempR_LO = (Type*)alloca(sizeof(Type)*nterms*3);

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;

  //do dumb averaging for q vector... this affects the velocity terms internally
  //TODO: use h.o. extrapolations here
  for(i = 0; i < neqn; i++){
    Qavg[i] = (qL[i] + qR[i])/2.0;
  }
  eqnset->ComputeAuxiliaryVariables(Qavg);

  std::vector<Int> gradientsLoc;
  eqnset->GetGradientsLocation(gradientsLoc);
  if(gradientsLoc.size() != (UInt)nterms){
    Abort << "WARNING: gradients location vector does not match number of terms in gradient";
  }

  memcpy(QL, qL, sizeof(Type)*nvars);
  memcpy(QR, qR, sizeof(Type)*nvars);
  eqnset->NativeToExtrapolated(QL);
  eqnset->NativeToExtrapolated(QR);

  //get turbulent viscosity
  Type tmut = 0.5*(mut[left_cv] + mut[right_cv]);

  //use directional derivatives to get gradient at face
  //TODO: generalize in case we want to use a non-midpoint edge CV
  for(i = 0; i < nterms*3; i++){
    grad[i] = 0.5*(gradL[i] + gradR[i]);
  }

  Type s2 = 0.0;
  Type qdots, dq;
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }
  for(j = 0; j < nterms; j++){
    qdots = DotProduct(dx, &grad[j*3]);
    Int loc = gradientsLoc[j];
    dq = (QR[loc] - QL[loc] - qdots)/s2;
    for(i = 0; i < 3; i++){
      grad[j*3 + i] += dq*dx[i];
    }
  }

  //call viscous flux routine
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempR);

  //call low order viscous flux routine, w/o extrapolated gradients
  for(i = 0; i < nterms*3; i++){
    grad[i] = 0.5*(gradL[i] + gradR[i]);
  }
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempR_LO);

  //compute difference in higher order and lower order fluxes
  for(i = 0; i < neqn; i++){
    tempR[i] -= tempR_LO[i];
  }

  Type* ETE_src = space->GetField("ETEVisc", FIELDS::STATE_NONE);
  
  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &ETE_src[left_cv*neqn];
  *ptrR = &ETE_src[right_cv*neqn];
  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempR[i];
  }

  return;
}

template <class Type>
void Bkernel_Viscous_Src(B_KERNEL_ARGS)
{
  Int i, j;
  EqnSet<Type>* eqnset = space->eqnset;
  Mesh<Type>* m = space->m;
  Type* mut = (Type*)custom;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nterms = space->grad->neqn;
  Type* qL = &space->q[left_cv*nvars];
  Type* qR = &space->q[right_cv*nvars];
  Type* QL = (Type*)alloca(sizeof(Type)*nvars);
  Type* QR = (Type*)alloca(sizeof(Type)*nvars);
  Type* Qavg = (Type*)alloca(sizeof(Type)*nvars);
  Type* xL;
  Type* xR;
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* gradL = space->qgrad + 3*nterms*left_cv;
  Type* grad = (Type*)alloca(sizeof(Type)*nterms*3);
  Type* tempL_LO = (Type*)alloca(sizeof(Type)*neqn);

  xL = m->cg + 3*left_cv;
  xR = m->cg + 3*right_cv;

  BoundaryConditions<Real>* bc = space->bc;
  Int bcNum = bc->bc_map[factag];
  Int bcId; 
  if(factag == 0){
    bcId = ParallelBoundary;
  }
  else{
    bcId = bc->bc_applied[bcNum];
  }  

  //use directional derivatives to get gradient at face
  //TODO: generalize in case we want to use a non-midpoint edge CV
  Type s2 = 0.0;
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }

  //do dumb averaging for q vector... this affects the velocity terms internally
  //TODO: use h.o. extrapolations here
  for(i = 0; i < neqn; i++){
    Qavg[i] = 0.5*(qL[i] + qR[i]);
  }
  eqnset->ComputeAuxiliaryVariables(Qavg);

  std::vector<Int> gradientsLoc;
  eqnset->GetGradientsLocation(gradientsLoc);
  if(gradientsLoc.size() != (UInt)nterms){
    Abort << "WARNING: gradients location vector does not match number of terms in gradient";
  }

  //get turbulent viscosity
  Type tmut;
  if(m->IsGhostNode(right_cv)){
    memcpy(QL, qL, sizeof(Type)*nvars);
    memcpy(QR, qR, sizeof(Type)*nvars);
    eqnset->NativeToExtrapolated(QL);
    eqnset->NativeToExtrapolated(QR);

    Type* gradR = space->qgrad + 3*nterms*right_cv;
    tmut = (mut[left_cv] + mut[right_cv])/2.0;
    for(i = 0; i < nterms*3; i++){
      grad[i] = 0.5*(gradL[i] + gradR[i]);
    }
    Type qdots, dq;
    for(j = 0; j < nterms; j++){
      qdots = DotProduct(dx, &grad[j*3]);
      Int loc = gradientsLoc[j];
      dq = (QR[loc] - QL[loc] - qdots)/s2;
      for(i = 0; i < 3; i++){
	grad[j*3 + i] += dq*dx[i];
      }
    }
  }
  else{
    //if we are on a boundary condition, this needs to be handled via a directional
    //derivative type rule, since our gradients do not include the boundary points
    tmut = mut[left_cv];
    //DirectionalDerivativeGrad(dx, s2, avec, qL, qR, grad, list, nterms);
    memcpy(grad, gradL, sizeof(Type)*3*nterms);
  }

  //call viscous flux routine
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempL);

  if(m->IsGhostNode(right_cv)){
    Type* gradR = space->qgrad + 3*nterms*right_cv;
    for(i = 0; i < nterms*3; i++){
      grad[i] = 0.5*(gradL[i] + gradR[i]);
    }
  }
  else{
    memcpy(grad, gradL, sizeof(Type)*3*nterms);
  }

  //call viscous flux routine with lower order flux
  eqnset->ViscousFlux(Qavg, grad, avec, tmut, tempL_LO);

  //compute difference in higher and lower order fluxes
  for(i = 0; i < neqn; i++){
    tempL[i] -= tempL_LO[i];
  }

  Type* ETE_src = space->GetField("ETEVisc", FIELDS::STATE_NONE);
  
  //set data necessary for driver scatter
  *size = neqn;
  *ptrL = &ETE_src[left_cv*neqn];
  //correct signs (tempspace)
  for(i = 0; i < neqn; i++){
    tempL[i] = -tempL[i];
  }

  return;
}
