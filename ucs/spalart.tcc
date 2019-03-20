#include "bc.h"
#include "eqnset.h"
#include "mesh.h"
#include "macros.h"
#include "dataInfo.h"
#include "solutionSpace.h"
#include "solutionField.h"
#include "bc_defines.h"
#include "powerLaw.h"

template <class Type>
Spalart<Type>::Spalart()
{
  //do not call this!!!
}

template <class Type>
Spalart<Type>::~Spalart()
{
  //nothing to see here, move along
}

template <class Type>
Spalart<Type>::Spalart(SolutionSpace<Type>* space) :
  TurbulenceModel<Type>(space)
{
  Mesh<Type>* m = space->m;
  this->neqn = 1;
  //initialize model constants
  sigma = 2.0/3.0;
  cb1 = 0.1355;
  cb2 = 0.622;
  kappa = 0.41;
  cw2 = 0.3;
  cw3 = 2.0;
  cv1 = 7.1;
  ct3 = 1.2;
  ct4 = 0.5;
  cw1 = cb1/(kappa*kappa) + (1.0 + cb2)/sigma;
  cw36 = cw3*cw3*cw3*cw3*cw3*cw3;
  cv13 = cv1*cv1*cv1;

  this->idata = new DataInfo(this->neqn, std::string("TurbulentVariables"));
  this->idata->AddScalar(0, "Turb-nu");
  this->idata->Verify();
  this->space->AddField(*this->idata, FIELDS::STATE_TIME, FIELDS::VAR_EVERYWHERE);
  this->tvar = this->space->GetFieldData("TurbulentVariables", FIELDS::STATE_NP1);
  this->tvarold = this->space->GetFieldData("TurbulentVariables", FIELDS::STATE_N);
  this->tvaroldm1 = this->space->GetFieldData("TurbulentVariables", FIELDS::STATE_NM1);

  DataInfo tdata(3, std::string("Grad-TurbNu"));
  tdata.AddVector(0, "Grad-TurbNu");
  tdata.Verify();
  this->space->AddField(tdata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  this->tgrad = this->space->GetFieldData("Grad-TurbNu", FIELDS::STATE_NONE);

  //allocate array for storing infinity values
  this->tvarinf = new Type[this->neqn];
  //set infinity values for turbulence variables
  SetTinf();
}

template <class Type>
void Spalart<Type>::SetTinf()
{
 //get nu_ref, same as nu_inf
  Param<Type>* param = this->space->param;
  Type nu_ref = param->ref_viscosity/param->ref_density;
  
  //initialize inifinity values for BCs use
  //set the variable to farfield value of the turbulence variables
  //spalart's paper says this can be anywhere from 3-5 nuinf
  //this->tvarinf[0] = 5.0*nu_ref;
  
  //from MSU
  this->tvarinf[0] = 1.341946;

  std::cout << "SPALART: nu_inf = " << this->tvarinf[0] << std::endl;
}

template <class Type>
void Spalart<Type>::Initialize()
{
  std::cout << "SPALART: Initializing turbulence model" << std::endl;

  EqnSet<Type>* eqnset = this->space->eqnset;
  Mesh<Type>* m = this->space->m;
  Param<Type>* param = this->space->param;
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  
  Type* mut = this->space->GetFieldData("mut", FIELDS::STATE_NONE);
  Type* dist = this->space->GetFieldData("wallDistance", FIELDS::STATE_NONE);
  Type* q = this->space->q;

  Type Re = param->Re;
  if(!param->useRestart){
    for(Int i = 0; i < nnode+gnode; i++){
      //this is in Spalart's paper
      Type nu_ref = param->ref_viscosity/param->ref_density;
      //this->tvar[i] = nu_ref / 10.0;
      
      this->tvar[i] = this->tvarinf[0];
      
      //NGC - can't seem to get the solver to run with such harsh initialization
      //      try using a linear variation from the spalart init to the freestream
      //if(real(dist[i]) <= 1.0){
      //  this->tvar[i] = (this->tvarinf[0] - nu_ref/10.0)/1.0*dist[i] + nu_ref/10.0;
      //}
      //else{
      //  this->tvar[i] = this->tvarinf[0];
      //}
      
      //sets up a power law distribution from zero to tvarinf
      //this->tvar[i] = PowerLawU(this->tvarinf[0], dist[i], Re);

      this->tvarold[i] = this->tvar[i];
      this->tvaroldm1[i] = this->tvar[i];
    }
  }

  //initialize flow field turbulent viscosity
  Int qnvars = eqnset->neqn + eqnset->nauxvars;
  for(Int i = 0; i < nnode+gnode; i++){
    Type* ql = &q[i*qnvars];
    Type rho = eqnset->GetDensity(ql);
    Type mu = eqnset->ComputeViscosity(ql);
    Type nu = mu/rho;
    mut[i] = ComputeEddyViscosity(rho, nu, i);
  }

  //allocate a crs object to take care of the solution duties
  this->crs.Init(nnode, gnode, this->neqn, m->ipsp, m->psp, this->space->p);
  this->limiter = new Limiter<Type>(this->space, this->tvar, this->tgrad, this->neqn, this->neqn, this->neqn, "spalart");
}

template <class Type>
void Spalart<Type>::BC_Kernel(B_KERNEL_ARGS)
{
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  BoundaryConditions<Real>* bc = space->bc;
  Mesh<Type>* m = space->m;
  Int bcType = bc->GetBCType(factag); 
  Type* tL = &turb->tvar[left_cv*turb->neqn];
  Type* tR = &turb->tvar[right_cv*turb->neqn];

  *ptrL = NULL;
  *size = 0;

  if(bcType == Proteus_ParallelBoundary){
    //do nothing this is a parallel updated edge/node
  }
  else if(bcType == Proteus_NoSlip){
    //this is a viscous surface, set the bc to zero
    *tR = 0.0;
    *tL = 0.0;
  }
  else if(bcType == Proteus_Symmetry || bcType == Proteus_ImpermeableWall){
    *tR = *tL;
  }
  else{
    //set the variable to farfield value of the turbulence variables
    //spalart's paper says this can be anywhere from 3-5 nuinf
    *tR = this->tvarinf[0];
  }
}

template <class Type>
void Spalart<Type>::BC_Jac_Kernel(B_KERNEL_ARGS)
{
  TurbulenceModel<Type>* turb = (TurbulenceModel<Type>*) custom;
  CRS<Type>* crs = (CRS<Type>*) &turb->crs;
  BoundaryConditions<Real>* bc = space->bc;
  Int neqn = turb->neqn;
  Int bcType = bc->GetBCType(factag); 

  *ptrL = NULL;
  *size = 0;

  if(bcType == Proteus_ParallelBoundary){
    //do nothing this is a parallel updated edge/node
  }
  else if(bcType == Proteus_NoSlip){
    //this is a viscous surface, set the bc to zero
    //we need to nuke the off-diagonal terms in the jacobian
    //as well as the residual since the bc is dirichlet on the wall
    crs->b[left_cv*neqn] = 0.0;
    crs->x[left_cv*neqn] = 0.0;
    crs->A->BlankRow(left_cv);
    Type* ptr = crs->A->GetPointer(left_cv, left_cv);
    //set diag to unity.. just to keep SGS/GMRES happy
    *ptr = 1.0;
  }
  else{
    //do nothing
  }
}

template <class Type>
void Spalart<Type>::Source(Type nu, Type d, Type* vgrad, Type* tvars, Type vol, 
			   Type* res, Type* jac)
{
  EqnSet<Type>* eqnset = this->space->eqnset;
  Type nut = *tvars;
  Type chi = nut/nu;
  Type chi2 = chi*chi;
  Type chi3 = chi2*chi;

  //trip function
  Type ft2 = ct3*exp(-ct4*chi2);
  Type d2 = d*d;
  Type prod, pi;
  Type dest, di;
  
  Type uy = vgrad[1];
  Type uz = vgrad[2];
  Type vx = vgrad[3];
  Type vz = vgrad[5];
  Type wx = vgrad[6];
  Type wy = vgrad[7];

  Type omega[3];
  omega[0] = wy - vz;
  omega[1] = uz - wx;
  omega[2] = vx - uy;

  Type Reinv = 1.0/eqnset->GetRe();
  Type Ret = eqnset->GetRe();

  Type fv1 = chi3/(chi3 + cv13);
  Type fv2 = 1.0 - chi/(1.0 + chi*fv1);
  Type s = Magnitude(omega) + (nut/(kappa*kappa*d2))*fv2*Reinv;

  //limiter for S
  //we limit S to be no smaller than 0.3*omega
  //this comes from the NASA Langley docs, private comm. with Spalart
  Type limitLow = 1.0e-12;
  s = MAX(MAX(s, 0.3*Magnitude(omega)), limitLow);

  //wake function
  Type r = MIN(10.0, nut/(s*kappa*kappa*d2)*Reinv);

  Type r6 = r*r*r;
  r6 = r6*r6;
  Type g = r + cw2*(r6 - r);
  Type g6 = g*g*g;
  g6 = g6*g6;
  Type fw = g*pow(((1.0 + cw36)/(g6 + cw36)), 1.0/6.0);

  //compute the production term
  pi = cb1*(1.0-ft2)*s;
  prod = pi*nut;

  //compute the destruction term
  di = (cw1*fw - cb1*ft2/(kappa*kappa))*Reinv*(nut/(d2));
  dest = di*nut;

  //compute ds/dnut, freeze fv2, I don't think the jacobians need the accuracy
  Type dsdnut = fv2*Reinv/(kappa*kappa*d2);
  //compute dP/dnut, freeze ft2
  Type dPdnut = pi*nut*dsdnut/s;
  //compute dD/dnut, freeze fw and ft2
  Type dDdnut = di;
  
  *res = (prod - dest)*vol;
  *jac = (MAX(0.0,-(pi-di)) + MAX(0.0, -(dPdnut - dDdnut)))*vol;
}

template <class Type>
void Spalart<Type>::Diffusive(Type nu, Type* tgrad, Type* tvarsL, Type* tvarsR, 
			      Type* avec, Type dgrad, Type* resL, Type* resR, 
			      Type* jacL, Type* jacR)
{
  Type nutL = tvarsL[0];
  Type nutR = tvarsR[0];
  Type nut = 0.5*(nutL + nutR);
  Type gdot = DotProduct(tgrad, avec);
  Type area = avec[3];
  
  Type Reinv = 1.0/this->space->eqnset->GetRe();

  //diffusive flux
  Type c1 = (1.0 + cb2)*nut + nu;
  //NOTE: for the second term -cb2.. we have assumed that
  //      nut is locally constant to deal with the grad(nut) * grad(nut) discretization
  *resL = c1*gdot - cb2*nutL*gdot;
  *resL *= 1.0/sigma*area*Reinv;

  *resR = c1*gdot - cb2*nutR*gdot;
  *resR *= 1.0/sigma*area*Reinv;
  
  //diffusive flux jacobians, computed using directional derivatives
  *jacL = c1*dgrad - cb2*nutL*dgrad;
  *jacL *= 1.0/sigma*area*Reinv;

  *jacR = c1*dgrad - cb2*nutR*dgrad;
  *jacR *= 1.0/sigma*area*Reinv;
}

template <class Type>
Type Spalart<Type>::ComputeEddyViscosity(Type rho, Type nu, Int node)
{
  Type nut = this->tvar[node*this->neqn];

  Type chi = nut/nu;

  //no turbulent viscosity
  if(nut == 0.0) return 0.0;

  Type chi3 = chi*chi*chi;
  Type fv1 = chi3 / (chi3 + cv13);
  Type mu_t = rho*nut*fv1;

  return mu_t;
}
