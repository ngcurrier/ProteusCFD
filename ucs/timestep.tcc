#include "macros.h"
#include "driver.h"
#include "eqnset.h"
#include "solutionSpace.h"

template <class Type>
Type ComputeTimesteps(SolutionSpace<Type>* space)
{
  Int i;
  Mesh<Type>* m = space->m;
  Param<Type>* param = space->param;
  EqnSet<Type>* eqnset = space->eqnset;
  Type dtmin = param->dt;
  Type CFL = param->GetCFL();
  Type VNN = param->VNN;
  Int nvars = eqnset->neqn+eqnset->nauxvars;
  Int nnode = m->GetNumNodes();

  Kernel<Type> Timestep(Kernel_Timestep);
  Kernel<Type> BTimestep(Bkernel_Timestep);

  Type* dt = space->GetField("timestep", FIELDS::STATE_NONE);
  // zero out the out timesteps
  for(i = 0; i < nnode; i++){
    dt[i] = 0.0;
  }
  // if we use local time stepping, compute the dt in the standard way
  if(param->useLocalTimeStepping){
    Driver(space, Timestep, nvars, (void*)dt);
    Bdriver(space, BTimestep, nvars, (void*)dt);
    dt[0] = CFL*(m->vol[0]/dt[0]);
    dtmin = dt[0];
    for(i = 1; i < nnode; i++){
      dt[i] = CFL*(m->vol[i]/dt[i]);
      //this line is the Von Neumann number diffusion stability limit
      //from D. Unrau "Viscous Airfoil Computations Using Local Preconditioning" 
      //Univ. of Toronto AIAA paper 2027 - yr. 1997
      if(param->enableVNN){
	dt[i] = MIN(dt[i], VNN*pow(m->vol[i], 2.0/3.0));
      }
      dtmin = MIN(dtmin, dt[i]);
    }
  }
  // if we don't use local time stepping nor pseudotime 
  // and the time step is specified, use that value everywhere
  else if(!param->pseudotimestepping && real(param->dt) > 0.0){
    for(i = 0; i < nnode; i++){
      dt[i] = param->dt;
    }
  }
  // if we don't use local time stepping and we use pseudotime
  // and the pseudo time step is specified, use that value everywhere
  else if(param->pseudotimestepping && real(param->dtau) > 0.0){
    for(i = 0; i < nnode; i++){
      dt[i] = param->dtau;
    }
  }
  // if we don't use local time stepping and the time step is 
  // NOT specified in pseudo nor real time, use the minimum allowable
  // timestep by the CFL condition in everycell
  // This keeps all the cells at the same time level but allows us
  // to march at the maximum allowable timestep by the CFL condition
  else{
    Driver(space, Timestep, nvars, (void*)dt);
    Bdriver(space, BTimestep, nvars, (void*)dt);
    dtmin = CFL*(m->vol[0]/dt[0]);
    for(i = 1; i < nnode; i++){
      dtmin = MIN(dtmin, CFL*(m->vol[i]/dt[i]));
      if(param->enableVNN){
	//same as above this is the von neumann number condition
	dtmin = MIN(dtmin, VNN*pow(m->vol[i], 2.0/3.0));
      }
    }
    for(i = 0; i < nnode; i++){
      dt[i] = dtmin;
    }
  }
  
  return (dtmin);
}

template <class Type>
void Kernel_Timestep(KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Type* dt = (Type*)custom;
  Type* beta = space->GetField("beta", FIELDS::STATE_NONE);
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type area = avec[3];
  Type maxeig;
  Type gamma = eqnset->param->gamma;
  //use passed in tempspace to compute average
  //of variables across edge
  for(i = 0; i < neqn; i++){
    tempL[i] = 0.5*(QL[i] + QR[i]); 
  }
  eqnset->ComputeAuxiliaryVariables(tempL);
  
  Type betaL = beta[left_cv];
  Type betaR = beta[right_cv];
  Type avbeta = 0.5*(betaL + betaR);

  maxeig = eqnset->MaxEigenvalue(tempL, avec, vdotn, gamma, avbeta); 

  //set data necessary for driver scatter
  *size = 1;
  *ptrL = &dt[left_cv];
  *ptrR = &dt[right_cv];
  tempL[0] = maxeig*area;
  tempR[0] = maxeig*area;

  return;
}

template <class Type>
void Bkernel_Timestep(B_KERNEL_ARGS)
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Type* dt = (Type*)custom;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* QL = &space->q[left_cv*nvars];
  Type* QR = &space->q[right_cv*nvars];
  Type area = avec[3];
  Type maxeig;
  Type gamma = eqnset->param->gamma;
  //use passed in tempspace to compute average
  //of variables across edge
  for(i = 0; i < neqn; i++){
    tempL[i] = 0.5*(QL[i] + QR[i]); 
  }
  eqnset->ComputeAuxiliaryVariables(tempL);

  Type* beta = space->GetField("beta", FIELDS::STATE_NONE);
  Type betaL = beta[left_cv];

  maxeig = eqnset->MaxEigenvalue(tempL, avec, vdotn, gamma, betaL); 

  //set data necessary for driver scatter
  *size = 1;
  *ptrL = &dt[left_cv];
  tempL[0] = maxeig*area;
  
  return;
}
