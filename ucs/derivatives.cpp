#include "derivatives.h"
#include "portFileio.h"
#include "move.h"
#include "solutionSpace.h"
#include "param.h"
#include "bc.h"
#include "parallel.h"
#include "macros.h"
#include "forces.h"
#include "mesh.h"  
#include "eqnset.h"
#include "param.h"
#include "create_functions.h"
#include "solve.h"
#include "solutionOperations.h"
#include "residual.h"
#include "mem_util.h"
#include "kernel.h"
#include "driver.h"
#include <iostream>

void Compute_dRdBeta_CTSE(Real* dRdB, SolutionSpace<Real>& space, Int beta)
{
  RCmplx perturb (0.0, 1.0e-11);

  Mesh<Real>* rm = space.m;

  Int nnode = rm->GetNumNodes();
  Int gnode = rm->GetNumParallelNodes();
  Int nbnode = rm->GetNumBoundaryNodes();

  std::ifstream fin;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;
  Param<RCmplx>* param = cspace.param;
  
  Int cnnode = cm->GetNumNodes();
  Int cgnode = cm->GetNumParallelNodes();
  Int cnbnode = cm->GetNumBoundaryNodes();

  RCmplx* q = cspace.q;

  Int neqn = ceqnset->neqn;
  Int nvars = neqn + ceqnset->nauxvars;

  Real* dxdb = new Real[nnode*3];

  Get_dXdB(space.param->path + space.param->spacename, dxdb, rm, beta);

  //perturb complex part by dxdb*perturb
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      cm->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  cp->UpdateXYZ(cm->xyz);
  //calculate metrics with new grid positions
  cm->CalcMetrics();
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&cspace);

  //perturb any parameters which we are getting sensitivities for
  Perturb_Param(beta, *cspace.param);
  //reset any interesting field which might depend on perturbations
  cspace.RefreshForParam();
  
  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(cspace.param->gaussianSource){
    cspace.gaussian->ApplyToResidual();
  }

  //update boundary conditions
  cspace.p->UpdateGeneralVectors(q, nvars);
  UpdateBCs(&cspace);
  cspace.p->UpdateGeneralVectors(q, nvars);

  //now, compute gradients and limiters
  if(param->sorder > 1 || param->viscous){
    cspace.grad->Compute();
    cspace.limiter->Compute(&cspace);
  }

  //compute spatial residual
  SpatialResidual(&cspace);
  ExtraSourceResidual(&cspace);
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < ceqnset->neqn; j++){
      dRdB[i*ceqnset->neqn + j] = imag(cspace.crs->b[i*ceqnset->neqn + j])/imag(perturb);
    }
  } 

  delete [] dxdb;
}

void Compute_dRdBeta_FD(Real* dRdB, SolutionSpace<Real>& space, Int beta)
{
  Real perturb = 1.0e-8;

  EqnSet<Real>* eqnset = space.eqnset;
  Mesh<Real>* m = space.m;
  PObj<Real>* p = space.p;
  Param<Real>* param = space.param;
  Real* q = space.q;
  Int i;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;

  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int nbnode = m->GetNumBoundaryNodes();

  Real* Qorig = new Real[nnode*nvars];
  Real* resOrig = new Real[nnode*neqn];
  Real* resUp = new Real[nnode*neqn];
  Real* resDown = new Real[nnode*neqn];

  //store original q
  for(i = 0; i < (nnode*nvars); i++){
    Qorig[i] = space.q[i];
  }

  UpdateBCs(&space);
  //compute spatial residual
  SpatialResidual(&space);
  ExtraSourceResidual(&space);
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < eqnset->neqn; j++){
      resOrig[i*eqnset->neqn + j] = space.crs->b[i*eqnset->neqn + j];
    }
  } 

  Real* dxdb = new Real[nnode*3];

  Get_dXdB(space.param->path+space.param->spacename, dxdb, m, beta);

  //perturb complex part by dxdb*perturb UP
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  //calculate metrics with new grid positions
  m->CalcMetrics();
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&space);

  //perturb parameters as required
  Perturb_Param(beta, *(space.param), 1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();

  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(space.param->gaussianSource){
    space.gaussian->ApplyToResidual();
  }

  //update boundary conditions
  UpdateBCs(&space);
  //now, compute gradients and limiters
  if(param->sorder > 1 || param->viscous){
    space.grad->Compute();
    space.limiter->Compute(&space);
  }

  //compute spatial residual
  SpatialResidual(&space);
  ExtraSourceResidual(&space);
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < eqnset->neqn; j++){
      resUp[i*eqnset->neqn + j] = (space.crs->b[i*eqnset->neqn + j]);
    }
  } 

  //perturb parameters as required, call twice (move back to original value then perturb down)
  Perturb_Param(beta, *(space.param), -1);
  Perturb_Param(beta, *(space.param), -1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();

  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(space.param->gaussianSource){
    space.gaussian->ApplyToResidual();
  }

  //perturb complex part by dxdb*perturb DOWN
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] -= 2.0*dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  //calculate metrics with new grid positions
  m->CalcMetrics();
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&space);
  //update boundary conditions
  UpdateBCs(&space);
  //now, compute gradients and limiters
  if(param->sorder > 1 || param->viscous){
    space.grad->Compute();
    space.limiter->Compute(&space);
  }
  //compute spatial residual
  SpatialResidual(&space);
  ExtraSourceResidual(&space);
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < eqnset->neqn; j++){
      resDown[i*eqnset->neqn + j] = (space.crs->b[i*eqnset->neqn + j]);
    }
  } 

  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < eqnset->neqn; j++){
      //up
      resUp[i*eqnset->neqn + j] -= resOrig[i*eqnset->neqn + j];
      resUp[i*eqnset->neqn + j] /= perturb;
      //down
      resDown[i*eqnset->neqn + j] -= resOrig[i*eqnset->neqn + j];
      resDown[i*eqnset->neqn + j] /= (-perturb);
      //central
      dRdB[i*eqnset->neqn + j] = (resUp[i*eqnset->neqn + j] + resDown[i*eqnset->neqn +j]) / 2.0;
    }
  }
  
  //reset metrics in case we need them
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  //recalculate metrics
  m->CalcMetrics();
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&space);
  //reset Qorig in case we need it
  for(i = 0; i < nnode*nvars; i++){
    space.q[i] = Qorig[i];
  }
 

  //go back to original values
  Perturb_Param(beta, *(space.param), 1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();

 //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(space.param->gaussianSource){
    space.gaussian->ApplyToResidual();
  }

 //update boundary conditions
  UpdateBCs(&space);
  //now, compute gradients and limiters
  if(param->sorder > 1 || param->viscous){
    space.grad->Compute();
    space.limiter->Compute(&space);
  }

  delete [] dxdb;

  delete [] resOrig;
  delete [] resUp;
  delete [] resDown;

  return;
}


void Compute_dQdBeta(Real* dQdB, SolutionSpace<Real>& space, Int beta)
{
  Int ndv;
  Real deltaDq;
  EqnSet<Real>* eqnset = space.eqnset;
  BoundaryConditions<Real>* bc = space.bc;
  Param<Real>* param = space.param;
  Mesh<Real>* m = space.m;
  Int nnode = m->GetNumNodes();
  Int neqn = eqnset->neqn;
  
  std::cout << "\n\nCOMPUTING dQ/dBeta\n" << std::endl;

  std::string fullpath = space.param->path + space.param->spacename;
  ndv = GetNdvDesignFile(fullpath);

  //allocate enough memory for dRdB
  Real* dRdB = new Real[nnode*neqn];

  if(1){
    Compute_dRdBeta_CTSE(dRdB, space, beta);
  }
  else{
    Compute_dRdBeta_FD(dRdB, space, beta);
  } 
 
  Compute_dRdQ(*space.crs, &space, 2);

  //copy drdb into res array in mesh
  memcpy(space.crs->b, dRdB, (nnode*neqn)*sizeof(Real));
  //set initial guess to zero
  space.crs->BlankX();

  if(param->designSolver == 0){
    //compute LU decomposition for the diagonal jacobians
    //SGS solver expects this
    space.crs->A->PrepareSGS();
    deltaDq = space.crs->SGS(param->designNsgs, NULL, NULL, NULL, true);
  }
  else if (param->designSolver == 1){
    deltaDq = space.crs->GMRES(param->designRestarts, param->designSearchDir, param->designPrecond, NULL, NULL, NULL);
  }  
  else{
    Abort << "Did not recognize solver type";
  }

  //copy dqdb out of dq array
  memcpy(dQdB, space.crs->x, (nnode*neqn)*sizeof(Real));
  std::cout.setf(std::ios::scientific);
  std::cout.precision(16);
  std::cout << "dQdBeta for beta " << beta << " D||linsolve-dq||: " 
	    << deltaDq << std::endl;
  
  delete [] dRdB;
  
  return;
}

void Compute_dQdBeta_CTSE(Real* dQdB, SolutionSpace<Real>& space, Int beta)
{
  RCmplx perturb(0.0, 1.0e-11);

  Mesh<Real>* rm = space.m;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;

  std::vector<SolutionSpaceBase<RCmplx>*> cSolSpaces;
  cSolSpaces.push_back(&cspace);

  SolutionOrdering<RCmplx> operations;
  std::string name = cspace.name;
  operations.Insert("Iterate " + name);
  operations.Finalize(cSolSpaces);

  Int cnnode = cm->GetNumNodes();

  Real* dxdb = new Real[cnnode*3];

  Get_dXdB(space.param->path+space.param->spacename, dxdb, rm, beta);

  //perturb complex part by dxdb*perturb
  for(Int i = 0; i < cnnode; i++){
    for(Int j = 0; j < 3; j++){
      cm->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  cp->UpdateXYZ(cm->xyz);
  //calculate metrics with new grid positions
  cm->CalcMetrics();
  
  //create a temporary operations vector with just a single space iterator on it
  //run solver
  Solve(cSolSpaces, operations);

  for(Int i = 0; i < cnnode; i++){
    for(Int j = 0; j < ceqnset->neqn; j++){
      dQdB[i*ceqnset->neqn + j] = 
	imag(cspace.q[i*(ceqnset->neqn+ceqnset->nauxvars) + j]) / imag(perturb);
    }
  }
  
  delete [] dxdb;

  return;
}

void Compute_dQdBeta_FD(Real* dQdB, SolutionSpace<Real>& space, Int beta)
{
  Real perturb = 1.0e-8;
  
  Int i ,j;
  EqnSet<Real>* eqnset = space.eqnset;
  Param<Real>* param = space.param;
  Mesh<Real>* m = space.m;
  PObj<Real>* p = space.p;

  Int nnode = m->GetNumNodes();

  std::vector<SolutionSpaceBase<Real>*> solSpaces;
  solSpaces.push_back(&space);

  //create a temporary operations vector with just a single space iterator on it
  SolutionOrdering<Real> operations;
  std::string name = space.name;
  operations.Insert("Iterate " + name);
  operations.Finalize(solSpaces);

  //compute dq/db with finite difference
  Real* dQdB_up = new Real[nnode*eqnset->neqn];
  Real* dQdB_down = new Real[nnode*eqnset->neqn];
  Real* Qorig = new Real[nnode*eqnset->neqn];

  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      Qorig[i*eqnset->neqn + j] = space.q[i*(eqnset->neqn+eqnset->nauxvars) + j];
    }
  }

  Real* dxdb = new Real[nnode*3];

  Get_dXdB(space.param->path+space.param->spacename, dxdb, m, beta);

  //perturb by dxdb*perturb UP
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  //recalculate metrics
  m->CalcMetrics();

  //perturb any parameters as required for sensitivities
  Perturb_Param(beta, *param, 1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();

  Solve(solSpaces, operations);
  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      dQdB_up[i*eqnset->neqn + j] = space.q[i*(eqnset->neqn+eqnset->nauxvars) + j];
    }
  }

  //perturb by dxdb*perturb DOWN
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] -= 2.0*dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  //recalculate metrics
  m->CalcMetrics();

  //perturb any parameters as required for sensitivities
  //back to original and then down
  Perturb_Param(beta, *param, -1);
  Perturb_Param(beta, *param, -1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();

  Solve(solSpaces, operations);
  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      dQdB_down[i*eqnset->neqn + j] = space.q[i*(eqnset->neqn+eqnset->nauxvars) + j];
    }
  }

  //compute differences
  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      //up
      dQdB_up[i*eqnset->neqn + j] -= Qorig[i*eqnset->neqn + j];
      dQdB_up[i*eqnset->neqn + j] /= perturb;
      //down
      dQdB_down[i*eqnset->neqn + j] -= Qorig[i*eqnset->neqn + j];
      dQdB_down[i*eqnset->neqn + j] /= -perturb;
    }
  }

  //average for central difference
  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      dQdB[i*eqnset->neqn + j] = 
	(dQdB_up[i*eqnset->neqn + j] + dQdB_down[i*eqnset->neqn + j]) / 2.0;
    }
  }

  //reset metrics in case we need them
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      m->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  p->UpdateXYZ(m->xyz);
  m->CalcMetrics();

  //go back to original value
  Perturb_Param(beta, *param, 1);
  //reset any interesting field which might depend on perturbations
  space.RefreshForParam();
  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(space.param->gaussianSource){
    space.gaussian->ApplyToResidual();
  }

  //update boundary conditions
  UpdateBCs(&space);

  //reset Qorig in case we need it
  for(i = 0; i < nnode; i++){
    for(j = 0; j < eqnset->neqn; j++){
      space.q[i*(eqnset->neqn+eqnset->nauxvars) + j] = Qorig[i*eqnset->neqn + j];
      eqnset->ComputeAuxiliaryVariables(&space.q[i*(eqnset->neqn+eqnset->nauxvars)]);
    }
  }

  delete [] Qorig;
  delete [] dQdB_up;
  delete [] dQdB_down;
  delete [] dxdb;

  return;
}

void Compute_dObjdQ(Real* dIdQ, SolutionSpace<Real>& space)
{
  RCmplx perturb (0.0, 1.0e-11);
  Int i, j;

  EqnSet<Real>* eqnset = space.eqnset;
  Mesh<Real>* rm = space.m;
  Param<Real>* param = space.param;
  BoundaryConditions<Real>* bc = space.bc;

  Int node;

  Int nnode = rm->GetNumNodes();
  Int gnode = rm->GetNumParallelNodes();
  Int nbnode = rm->GetNumBoundaryNodes();
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  MPI_Datatype mpit = MPI_GetType(nnode);

  //get global number of nodes since processes need to be synced
  //build lookup maps from local to global id
  Int* globalToLocal;
  Int* localToGlobal;
  Int globalNodes = space.p->BuildGlobalMaps(&localToGlobal, &globalToLocal);

  RCmplx Obj;

  std::cout << "\n\nCOMPUTING dI/dQ\n\n" << std::endl;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;

  //init dIdQ
  MemBlank(dIdQ, neqn*nnode);

  //compute surface areas once to ensure they are valid
  ComputeSurfaceAreas(&cspace, 0);

  for(i = 0; i < globalNodes; i++){
    node = globalToLocal[i];
    for(j = 0; j < neqn; j++){
      //perturb q value if node is local to process
      if(node >= 0){
	cspace.q[node*nvars + j] += perturb;
	ceqnset->ComputeAuxiliaryVariables(&cspace.q[node*nvars + 0]);
      }
      //compute perturbed objective function
      //no need for a call to compute_surface_areas here since we aren't changing them
      cp->UpdateGeneralVectors(cspace.q, nvars);
      Obj = Compute_Obj_Function(cspace);
      if(node >= 0){
	if(node < nnode){
	  dIdQ[node*neqn + j] = imag(Obj) / imag(perturb);
	}
	//reset q value
	cspace.q[node*nvars + j] = space.q[node*nvars + j];
	//compute aux vars one last time for all values reset to original
	ceqnset->ComputeAuxiliaryVariables(&cspace.q[node*nvars + 0]);
      }
    }
  }

  delete [] localToGlobal;
  delete [] globalToLocal;

  return;
}

void Compute_Adjoint_Variables(Real* lambda, SolutionSpace<Real>& space)
{
  Real deltaDq; 
  Mesh<Real>* m = space.m;
  Param<Real>* param = space.param;
  EqnSet<Real>* eqnset = space.eqnset;

  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int neqn = eqnset->neqn;

  Real* dIdQ = new Real[nnode*neqn];
  Compute_dObjdQ(dIdQ, space);
 
  std::cout << "Adjoint solve - attempting to allocate dRdQ transpose" << std::endl;
  CRS<Real>* dRdQT = new CRS<Real>;
  dRdQT->Init(nnode, gnode, eqnset->neqn, m->ipsp, m->psp, m->p);

  Compute_dRdQ_Transpose(*dRdQT, &space, 2);

  //copy dIdQ into res array in mesh
  memcpy(space.crs->b, dIdQ, (nnode*neqn)*sizeof(Real));

  //set initial guess to zero
  dRdQT->BlankX();

  if(param->designSolver == 0){
    //compute LU decomposition for the diagonal jacobians
    //SGS solver expects this
    dRdQT->A->PrepareSGS();
    deltaDq = dRdQT->SGS(param->designNsgs, NULL, NULL, NULL, true);
  }
  else if (param->designSolver == 1){
    deltaDq = dRdQT->GMRES(param->designRestarts, param->designSearchDir, param->designPrecond, NULL, NULL, NULL);
  }  
  else{
    Abort << "Did not recognize solver type";
  }

  //copy adjoint variables out
  memcpy(lambda, space.crs->x, (nnode*neqn)*sizeof(Real));

  std::cout.setf(std::ios::scientific);
  std::cout.precision(16);
  std::cout << "Adjoint solve - D||linsolve-res||: "  << deltaDq << std::endl;
  
  delete [] dIdQ;

  delete dRdQT;

  return;
}

void Compute_Adjoint_Variables_II(Real* lambda, SolutionSpace<Real>& space)
{
  Real deltaDq; 
  Mesh<Real>* m = space.m;
  Param<Real>* param = space.param;
  EqnSet<Real>* eqnset = space.eqnset;

  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int neqn = eqnset->neqn;

  Real* dIdQ = new Real[nnode*neqn];
  Real* prod = new Real[(nnode+gnode)*neqn];
  Compute_dObjdQ(dIdQ, space);

  //we have to transpose our left hand side operator as well
  ComputeJacobians(&space);
  space.crs->A->CRSTranspose();
  space.crs->A->PrepareSGS();

  std::cout << "Adjoint solve - attempting to allocate dRdQ transpose" << std::endl;
  CRS<Real>* dRdQT = new CRS<Real>;
  dRdQT->Init(nnode, gnode, eqnset->neqn, m->ipsp, m->psp, m->p);

  Compute_dRdQ_Transpose(*dRdQT, &space, 2);
  MemSet(lambda, 0.0, (nnode+gnode)*neqn);

  for(Int iter = 0; iter < space.temporalControl.nSteps; iter++){

    //compute matrix vector product of dRdQ^T and current lambda vector
    dRdQT->MatVecMultiply(dRdQT->A, lambda, prod);

    for(Int i = 0; i < nnode*neqn; i++){
      space.crs->b[i] = dIdQ[i] - prod[i];
    }
    
    //compute the global residual for the linear system
    Real residGlobal = ParallelL2Norm(space.p, space.crs->b, nnode*neqn);
    space.residOutFile << iter << ": ||II_res||: " << residGlobal;

    //set initial guess to zero
    space.crs->BlankX();
    Real deltaDq = space.crs->SGS(param->designNsgs, NULL, NULL, NULL, false);
    space.residOutFile << " ||II_dq||: " << deltaDq << std::endl;

    std::cout << iter << ": ||II_res||: " << residGlobal << " ||II_dq||: " << deltaDq << std::endl;

    //update our estimate of lambda
    for(Int i = 0; i < nnode*neqn; i++){
      lambda[i] += space.crs->x[i];
    }
    
    //we have to do a parallel sync for the matrix vector multiply
    space.p->UpdateGeneralVectors(lambda, neqn);
  }

  delete [] dIdQ;
  delete [] prod;

  delete dRdQT;

  return;
}


Real Compute_dObjdBeta_Adjoint(Real* lambda, SolutionSpace<Real>& space, Int beta)
{
  Int i;
  Real dobjdbeta;
  MPI_Datatype mpit;
  Mesh<Real>* m = space.m;
  Param<Real>* param = space.param;
  BoundaryConditions<Real>* bc = space.bc;
  Int neqn = space.eqnset->neqn;
  
  Int nnode = m->GetNumNodes();

  Real* dRdB = new Real[nnode*neqn];

  Compute_dRdBeta_CTSE(dRdB, space, beta);
  
  dobjdbeta = 0.0;
  for(i = 0; i < nnode*neqn; i++){
    dobjdbeta += dRdB[i]*lambda[i];
  }

  //get mpi datatype to send
  mpit = MPI_GetType(dobjdbeta);
  //sum up across all processes
  MPI_Allreduce(MPI_IN_PLACE, &dobjdbeta, 1, mpit, MPI_SUM, MPI_COMM_WORLD);

  //add dI/dx*dx/dB contributions
  dobjdbeta += Compute_dObjdX_dXdBeta(space, beta);

  delete [] dRdB;
 
  return dobjdbeta;
}


Real Compute_dObjdQ_dQdBeta(Real* dQdB, SolutionSpace<Real>& space)
{
  RCmplx perturb (0.0, 1.0e-11);
  RCmplx Obj;

  EqnSet<Real>* eqnset = space.eqnset;
  Mesh<Real>* m = space.m;
  Param<Real>* param = space.param;
  BoundaryConditions<Real>* bc = space.bc;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;

  Int cnnode = cm->GetNumNodes();

  Real dObjdBeta;
  Int i, j;
  Int neqn = eqnset->neqn;
  Int nauxvars = eqnset->nauxvars;

  //replace all Q values in mesh by dQdB*perturb and compute cost function
  for(i = 0; i < cnnode; i++){
    for(j = 0; j < neqn; j++){
      cspace.q[i*(neqn+nauxvars) + j] += dQdB[i*neqn + j]*perturb;
    }
    ceqnset->ComputeAuxiliaryVariables(&cspace.q[i*(neqn+nauxvars)]);
  }
  cp->UpdateGeneralVectors(cspace.q, neqn+nauxvars);
  UpdateBCs(&cspace);
  cp->UpdateGeneralVectors(cspace.q, neqn+nauxvars);

  ComputeSurfaceAreas(&cspace, 0);
  Obj = Compute_Obj_Function(cspace);
  
  dObjdBeta = imag(Obj)/imag(perturb);

  return dObjdBeta;
}

Real Compute_dObjdX_dXdBeta(SolutionSpace<Real>& space, Int beta)
{
  RCmplx perturb (0.0, 1.0e-11);

  Real dObjdBeta;
  RCmplx Obj;

  Mesh<Real>* m = space.m;
  EqnSet<Real>* eqnset = space.eqnset;
  Param<Real>* param = space.param;
  BoundaryConditions<Real>* bc = space.bc;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;

  Int nnode = m->GetNumNodes();
  Int cnnode = cm->GetNumNodes();

  Real* dxdb = new Real[nnode*3];

  Get_dXdB(space.param->path+space.param->spacename, dxdb, m, beta);

  //perturb complex part by dxdb*perturb
  for(Int i = 0; i < cnnode; i++){
    for(Int j = 0; j < 3; j++){
      cm->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  cp->UpdateXYZ(cm->xyz);
  //calculate metrics with new grid positions
  cm->CalcMetrics();

  ComputeSurfaceAreas(&cspace, 0);
  Obj = Compute_Obj_Function(cspace);

  dObjdBeta = imag(Obj) / imag(perturb);
  
  delete [] dxdb;

  return dObjdBeta;
}

Real Compute_dObjdBeta_Direct(SolutionSpace<Real>& space, Int beta)
{
  Real dObjdBeta;
  Mesh<Real>* m = space.m;
  EqnSet<Real>* eqnset = space.eqnset;

  Int neqn = eqnset->neqn;
  Int nnode = m->GetNumNodes();

  Real* dQdB = new Real[nnode*neqn]; 
  
  if(1){
    //this uses the incremental iterative method with matrix free 
    //matrix vector product
    Compute_dQdBeta_II(dQdB, space, beta);
  }
  else if(0){
    //this solves a linear system that depends on the full linearization dR/dQ
    //to be exact for it to function
    Compute_dQdBeta(dQdB, space, beta);
  }
  else if(0){
    //this computes the solution in CTSE space
    Compute_dQdBeta_CTSE(dQdB, space, beta);
  }
  else{
    //this uses standard finite differences and requires two more solutions
    Compute_dQdBeta_FD(dQdB, space, beta);
  }
  dObjdBeta = 0.0;
  dObjdBeta += Compute_dObjdQ_dQdBeta(dQdB, space);
  dObjdBeta += Compute_dObjdX_dXdBeta(space, beta);

  delete [] dQdB;

  return dObjdBeta;
}

void Compute_dXdB(SolutionSpace<Real>& space)
{
  Mesh<Real>* m = space.m;
  BoundaryConditions<Real>* bc = space.bc;

  RCmplx h (0.0, 1.0e-11);

  Int i, j;

  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int nbnode = m->GetNumBoundaryNodes();

  Int beta;

  //this contains dxdb dydx and dzdb
  Real* dxdb = new Real[nnode*3];
  RCmplx* dxdbSurfC = new RCmplx[nnode*3];

  Mesh<RCmplx> cm(*m);
  PObj<RCmplx> cp;
  cm.SetParallelPointer(&cp);
  cp.BuildCommMaps(&cm);
  //check for parallel comm sanity
  cp.CheckSanityCoords(cm.xyz);

  Int cnnode = cm.GetNumNodes();
  
  std::cout << "\n\nCOMPUTING dX/dB\n" << std::endl;

  std::string fullpath = space.param->path + space.param->spacename;
  Int ndv = GetNdvDesignFile(fullpath);
  std::cout << "\n\nFOUND " << ndv << " design variables\n" << std::endl;

  std::ofstream fout;
  std::stringstream ss;
  ss.clear();
  ss.str("");
  ss << cp.GetRank();
  std::string dxdbFilename = space.param->path+space.param->spacename + "-DxDb." + (ss.str()) + ".dat";
  fout.open(dxdbFilename.c_str());
  fout << std::setprecision(16);

  Int np = m->p->GetNp();
  Int rank = m->p->GetRank();
  fout << np << " " << rank << " " << cnnode << std::endl;

  for(beta = 0; beta < ndv; beta++){
    std::cout << "PERTURBING BETA: " << beta << std::endl;
    Compute_dXdB_Surface(space.param->path + space.param->spacename, m, bc, dxdb, beta);
    //perturb points by h * dxdb_surface
    for(i = 0; i < nnode; i++){
      for(j = 0; j < 3; j++){
	dxdbSurfC[i*3 + j] = dxdb[i*3 + j]*h;
      }
    }      
    Int smoothingPasses = 1000;
    MoveMeshLinearElastic(&cm, bc, dxdbSurfC, smoothingPasses);
  
    //compute dxdb
    for(i = 0; i < cnnode*3; i++){
      dxdb[i] = imag(cm.xyz[i])/imag(h);
    }
    //write dxdb dxdb dzdb to file
    for(i = 0; i < cnnode; i++){
      for(j = 0; j < 3; j++){
	fout << dxdb[i*3 + j] << " " ;
      }
      fout << std::endl;
    }
    //reset xyz coords for next pass
    for(i = 0; i < nnode*3; i++){
      cm.xyz[i] = m->xyz[i];
    }
    cp.UpdateXYZ(cm.xyz);

  }

  fout.close();

  delete [] dxdb;
  delete [] dxdbSurfC;

  return;
}

void Get_dXdB(std::string path, Real* dxdb, Mesh<Real>* m, Int beta)
{
  std::ifstream fin;
  std::stringstream ss;

  Int i, b;

  PObj<Real>* p = m->p;
  Int rank = p->GetRank();
  Int np = p->GetNp();

  Int rankRead, npRead, nnodeRead;
  Int nnode = m->GetNumNodes();

  ss.clear();
  ss << rank;

  //read points to move from file and dxdb for them
  std::string Filename = path + "-DxDb." + (ss.str()) + ".dat";
  fin.open(Filename.c_str());
  if(!fin.is_open()){
    std::cerr << "WARNING: opening file failed --> " << Filename << std::endl;
    return;
  }
  
  //read in header to check for mismatched process numbers
  fin >> npRead;
  fin >> rankRead;
  fin >> nnodeRead;

  //check that the state is sane
  if(npRead != np){
    std::cerr << "WARNING: Get_dXdB() Mesh movement sensitivities were run with a " 
	      << "nonmatching number of processes, files not synced!" << std::endl;
  }
  if(rankRead != rank){
    std::cerr << "WARNING: Get_dXdB() Mesh movement sensitivities opened with wrong rank!"
       << std::endl;
  }
  if(nnodeRead != nnode){
    std::cerr << "WARNING: Get_dXdB() Mesh movement sensitivities opened with wrong number "
	      << "of nodes!" << std::endl;
  }


  //read in globalnodes for each beta up until we have read
  //the beta we are looking for
  for(b = 0; b <= beta; b++){
    for(i = 0; i < nnode*3; i++){
      fin >> dxdb[i];
    }
  }

  fin.close();

  return;
}

Real Compute_dObjdBeta_CTSE(SolutionSpace<Real>& space, SolutionOrdering<Real>& operations, Int beta)
{
  Int i;  
  Real h = 1.0e-11;
  RCmplx I(0.0, 1.0);
  RCmplx perturb = h*I;

  Real dObjdBeta;
  RCmplx Obj;
  Mesh<Real>* rm = space.m;
  Param<Real>* param = space.param;
  BoundaryConditions<Real>* bc = space.bc;
  PObj<Real>* p = space.p;
  Int nnode = rm->GetNumNodes();

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;
  Int cnnode = cm->GetNumNodes();

  std::cout << "\n\nCOMPUTING dObj/dBeta_CTSE\n" << std::endl;

  Real* dxSurf = new Real[nnode*3];
  RCmplx* dxSurfC = new RCmplx[nnode*3];
  Compute_dXdB_Surface(space.param->path + space.param->spacename, rm, bc, dxSurf, beta);

  std::vector<SolutionSpaceBase<RCmplx>*> cSpaces;
  cSpaces.push_back(&cspace);

  for(i = 0; i < nnode*3; i++){
    //perturb complex part by displacement
    dxSurfC[i] = dxSurf[i]*perturb;
  }
  Int smoothingPasses = 1000;
  MoveMeshLinearElastic(cm, bc, dxSurfC, smoothingPasses);

  cm->CalcMetrics();

  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&cspace);

  SolutionOrdering<RCmplx> cOperations;
  cOperations = operations;
  cOperations.Finalize(cSpaces);

  //perturb parameters as required
  Perturb_Param(beta, *cspace.param);
  //reset any interesting field which might depend on perturbations
  cspace.RefreshForParam();

  //run solver
  Solve(cSpaces, cOperations);

  Obj = Compute_Obj_Function(cspace);
  dObjdBeta = imag(Obj) / h;
  
  std::cout << "\n\nComplex objective function value: " << Obj << std::endl;

  delete [] dxSurf;
  delete [] dxSurfC;

  return dObjdBeta;
}


void Compute_dRdQ_Product_MatrixFree(SolutionSpace<RCmplx>& cspace, Real* vector, Real* prod)
{
  RCmplx perturb(0.0, 1.0e-11);
  Mesh<RCmplx>* cm = cspace.m;
  Param<RCmplx>* param = cspace.param;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;
  RCmplx* q = cspace.q;
  PObj<RCmplx>* p = cspace.p;

  Int cnnode = cm->GetNumNodes();
  Int cgnode = cm->GetNumParallelNodes();
  Int cnbnode = cm->GetNumBoundaryNodes();

  Int neqn = ceqnset->neqn;
  Int nvars = ceqnset->nauxvars + neqn;

  for(Int i = 0; i < cnnode; i++){
    RCmplx* iq = &q[i*nvars + 0];
    for(Int j = 0; j < neqn; j++){
      iq[j] += vector[i*neqn + j]*perturb;
    }
    ceqnset->ComputeAuxiliaryVariables(iq);
  }
  p->UpdateGeneralVectors(q, nvars);

  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(cspace.param->gaussianSource){
    cspace.gaussian->ApplyToResidual();
  }

  UpdateBCs(&cspace);
  p->UpdateGeneralVectors(q, nvars);
  
  //now, compute gradients and limiters
  if(param->sorder > 1 || param->viscous){
    cspace.grad->Compute();
    cspace.limiter->Compute(&cspace);
  }

  //it is assumed that this routine is going to be called iteratively,
  //therefore, calling the turbulence model here, will slowly converge the
  //turbulence model towards the correct sensitivity
  //This should theoretically be converged here at every call to this routine,
  //but that is horribly costly.... so we cheat
  //WARNING: ---- for now, we assume frozen turbulence models
  if(param->viscous && false){
    cspace.turb->Compute();
    if(param->gcl){
      std::cout << "MATRIX FREE: COMPUTE GCL FOR TURB MODEL!" << std::endl;
    }
  }

  //Now compute residual, this only works for spatial residual right now
  SpatialResidual(&cspace);
  ExtraSourceResidual(&cspace);

  //Now do the CTSE derivative part
  for(Int i = 0; i < cnnode; i++){
    for(Int j = 0; j < neqn; j++){
      prod[i*neqn + j] = imag(cspace.crs->b[i*neqn + j])/imag(perturb);
    }
  }

  //re-zero the complex part of the q vector in case we reuse the cspace
  for(Int i = 0; i < cnnode; i++){
    RCmplx* iq = &q[i*nvars + 0];
    for(Int j = 0; j < neqn; j++){
      iq[j] = real(iq[j]);
    }
    ceqnset->ComputeAuxiliaryVariables(iq);
  }
  p->UpdateGeneralVectors(q, nvars);

  return;
}

void Compute_dQdBeta_II(Real* dQdB, SolutionSpace<Real>& space, Int beta)
{  
  Mesh<Real>* m = space.m;
  EqnSet<Real>* eqnset = space.eqnset;
  Param<Real>* param = space.param;
  Int neqn = eqnset->neqn;
  Int nnode = m->GetNumNodes();
  
  std::cout << "\n\nCOMPUTING dQ/dBeta - Incremental Iterative - Beta: " << beta << "\n";

  //allocate enough memory for dRdB
  Real* dRdB = new Real[nnode*neqn];
  Compute_dRdBeta_CTSE(dRdB, space, beta);

  SolutionSpace<RCmplx> cspace(space);

  //zero dQdB
  MemSet(dQdB, 0.0, nnode*neqn);
  
  //Get the direct computed dQdB on boundaries which are dirchlet(TODO) or viscous
  Int* vnodeslist;
  Int vnodes = GetNodesOnBCType(m, space.bc, &vnodeslist, NoSlip);
  Real* vnodesdQdB = new Real[vnodes*neqn];
  ComputedQdBeta_HardsetBC(vnodes, vnodeslist, vnodesdQdB, space, beta);

  //now hardset for the hard bcs
  for(Int i = 0; i < vnodes; i++){
    Int node = vnodeslist[i];
    for(Int j = 0; j < neqn; j++){
      //only non-zero entries should be hardset, the rest should float above
      if(Abs(vnodesdQdB[i*neqn + j]) >= 1.0e-16){
	dQdB[node*neqn + j] = vnodesdQdB[i*neqn + j];
      }
    }
  }

  for(Int it = 0; it < space.temporalControl.nSteps; it++){
    //build RHS of our problem
    //
    
    //compute the rhs contribution of dR/dQ*dQ/dB
    Compute_dRdQ_Product_MatrixFree(cspace, dQdB, space.crs->b);
    
    //flip the sign and add the dR/dB contribution
    for(Int i = 0; i < nnode*neqn; i++){
      space.crs->b[i] = -space.crs->b[i] - dRdB[i];
    }

    //this is a hook to modify any residuals directly due to boundary
    //conditions, most notably hardset viscous BCs (zeros appropriate entries)
    Kernel<Real> ResModifyBC(Bkernel_BC_Res_Modify);
    BdriverNoScatter(&space, ResModifyBC, eqnset->neqn, NULL);
    
    //compute the global residual for the linear system
    Real residGlobal = ParallelL2Norm(space.p, space.crs->b, nnode*neqn);
    space.residOutFile << it << ": ||II_res||: " << residGlobal;

    //set initial guess to zero
    space.crs->BlankX();
    Real deltaDq = space.crs->SGS(param->designNsgs, NULL, NULL, NULL, true);
    space.residOutFile << " ||II_sgs_res||: " << deltaDq << std::endl;

    std::cout << it << ": ||II_res||: " << residGlobal << " ||II_dq||: " << deltaDq << std::endl;

    //update our estimate of dQdB
    for(Int i = 0; i < nnode*neqn; i++){
      dQdB[i] -= space.crs->x[i];
    }

    //now hardset for the hard bcs
    for(Int i = 0; i < vnodes; i++){
      Int node = vnodeslist[i];
      for(Int j = 0; j < neqn; j++){
	//only non-zero entries should be hardset, the rest should float above
	if(Abs(vnodesdQdB[i*neqn + j]) >= 1.0e-16){
	  dQdB[node*neqn + j] = vnodesdQdB[i*neqn + j];
	}
      }
    }
  }

  delete [] dRdB;
  delete [] vnodeslist;
  delete [] vnodesdQdB;
}


void Perturb_Param(Int beta, Param<RCmplx>& param)
{ 
  Real h = 1.0e-11;
  RCmplx I(0.0, 1.0);
  RCmplx perturb = h*I;
  
  std::string designName = param.path+param.spacename;
  Int ndv = GetNdvDesignFile(designName);
  Real* x = new Real[ndv];
  Real* bounds = new Real[ndv*2];
  Real f;
  Real* grad = new Real[ndv];
  Int* dvType = new Int[ndv];

  ReadDesignFile(designName, &ndv, x, bounds, &f, grad, dvType); 

  Int parType = dvType[beta];

  if(parType == 0 || parType == 1 || parType == 2){
    //direct point movement x,y,z
  }
  else if(parType == 3 || parType == 4 || parType == 5){
    //hicks-henne functions
  }
  else if(parType == 6 || parType == 7 || parType == 8){
    //boundary movement x,y,z
  }
  else if(parType == 9){
    //x location gaussian plume
    param.gaussianXloc += perturb;
  }
  else if(parType == 10){
    //y location gaussian plume
    param.gaussianYloc += perturb;
  }
  else{
    //should never be here
    std::stringstream ss;
    ss << "In Perturb_Param() and design parameter type ";
    ss << parType;
    ss << " not found!";
    
    Abort << ss.str();
  }

  param.PostCompute();
  param.PrintSolverParams();
  
  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;
}

void Perturb_Param(Int beta, Param<Real>& param, Int sgn)
{ 
  Real perturb = 1.0e-8*(Real)sgn;

  std::string designName = param.path+param.spacename;
  Int ndv = GetNdvDesignFile(designName);
  Real* x = new Real[ndv];
  Real* bounds = new Real[ndv*2];
  Real f;
  Real* grad = new Real[ndv];
  Int* dvType = new Int[ndv];

  ReadDesignFile(designName, &ndv, x, bounds, &f, grad, dvType); 

  Int parType = dvType[beta];

  if(parType == 0 || parType == 1 || parType == 2){
    //direct point movement x,y,z
  }
  else if(parType == 3 || parType == 4 || parType == 5){
    //hicks-henne functions
  }
  else if(parType == 6 || parType == 7 || parType == 8){
    //boundary movement x,y,z
  }
  else if(parType == 9){
    //x location gaussian plume
    param.gaussianXloc += perturb;
  }
  else if(parType == 10){
    //y location gaussian plume
    param.gaussianYloc += perturb;
  }
  else{
    //should never be here
    std::stringstream ss;
    ss << "In Perturb_Param() and design parameter type ";
    ss << parType;
    ss << " not found!";
    
    Abort << ss.str();
  }

  param.PostCompute();
  param.PrintSolverParams();

  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;
}

void ComputedQdBeta_HardsetBC(Int nnodes, Int* nodes, Real* dQdB, const SolutionSpace<Real>& space, Int beta)
{
  Real h = 1.0e-11;
  RCmplx I(0.0, 1.0);
  RCmplx perturb = h*I;

  Mesh<Real>* rm = space.m;

  SolutionSpace<RCmplx> cspace(space);
  Mesh<RCmplx>* cm = cspace.m;
  PObj<RCmplx>* cp = cspace.p;
  EqnSet<RCmplx>* ceqnset = cspace.eqnset;
  Param<RCmplx>* param = cspace.param;
  RCmplx* q = cspace.q;

  Int nnode = rm->GetNumNodes();

  Int neqn = ceqnset->neqn;
  Int nvars = neqn + ceqnset->nauxvars;

  Real* dxdb = new Real[nnode*3];

  Get_dXdB(space.param->path + space.param->spacename, dxdb, rm, beta);

  //perturb complex part by dxdb*perturb
  for(Int i = 0; i < nnode; i++){
    for(Int j = 0; j < 3; j++){
      cm->xyz[i*3 + j] += dxdb[i*3 + j]*perturb;
    }
  }
  //update xyz coords
  cp->UpdateXYZ(cm->xyz);
  //calculate metrics with new grid positions
  cm->CalcMetrics();
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(&cspace);

  //perturb any parameters which we are getting sensitivities for
  Perturb_Param(beta, *cspace.param);
  //reset any interesting field which might depend on perturbations
  cspace.RefreshForParam();
  
  //the gaussian source is sometimes used to modify boundary velocities, etc.
  //pre-compute it for bc call
  if(cspace.param->gaussianSource){
    cspace.gaussian->ApplyToResidual();
  }

  //update boundary conditions
  cspace.p->UpdateGeneralVectors(q, nvars);
  UpdateBCs(&cspace);
  cspace.p->UpdateGeneralVectors(q, nvars);

  for(Int i = 0; i < nnodes; i++){
    Int node = nodes[i];
    for(Int j = 0; j < neqn; j++){
      dQdB[i*neqn + j] = imag(q[node*nvars + j])/h;
      std::cout << "dqdb: " << i << ": -" << j << " " <<  dQdB[i*neqn + j] << std::endl;
    }
  }

}
