#include "param.h"
#include "macros.h"
#include "matrix.h"
#include "solutionSpace.h"
#include "solutionField.h"

template <class Type>
HeatTransferEqnSet<Type>::HeatTransferEqnSet(SolutionSpace<Type>* space, Param<Type>* p)
{
  this->neqn = 1;
  this->nauxvars = 0;
  this->space = space;
  //variable set is conservative
  this->varsConservative = 1;

  //set data required
  this->idata = new DataInfo(this->neqn+this->nauxvars, "variableQ");
  this->idata->AddScalar(0, "Temperature");
  this->idata->Verify();

  //set gradients required
  this->gdata = new DataInfo((this->neqn)*3, "gradVariableQ");
  this->gdata->AddVector(0*3, "Grad-Temperature");
  this->gdata->Verify();
  
  this->Qinf = new Type[this->neqn + this->nauxvars];
  this->param = p;
}

template <class Type>
HeatTransferEqnSet<Type>::~HeatTransferEqnSet()
{
  delete [] this->Qinf;
  return;
}

template <class Type>
void HeatTransferEqnSet<Type>::InitEqnSet()
{

  this->UpdateQinf();

  //allocate solution memory
  this->space->AddField(*this->idata, FIELDS::STATE_TIME, FIELDS::VAR_EVERYWHERE);
  SolutionField<Type> & field = this->space->GetField("variableQ");
  this->space->q = field.GetData(FIELDS::STATE_NP1);
  this->space->qold = field.GetData(FIELDS::STATE_N);
  this->space->qoldm1 = field.GetData(FIELDS::STATE_NM1);

  //allocate solution memory for gradients we need
  this->param->viscous = true;
  this->space->AddField(*this->gdata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  SolutionField<Type> & gfield  = this->space->GetField("gradVariableQ");
  this->space->qgrad = gfield.GetData(FIELDS::STATE_NONE);
}


template <class Type>
void HeatTransferEqnSet<Type>::UpdateQinf()
{
  //assume everything starts at equilibrium with reference temperature
  this->Qinf[0] = 1.0;
}

template <class Type>
void HeatTransferEqnSet<Type>::SetInitialConditions()
{
   Int i, j;
   Int neqn = this->neqn;
   Int nauxvars = this->nauxvars;
   Mesh<Type>* m = this->space->m;
   Int nnode = m->GetNumNodes();
   Int gnode = m->GetNumParallelNodes();
   Int nbnode = m->GetNumBoundaryNodes();

   //make sure qinf is updated
   this->UpdateQinf();
   std::cout.setf(std::ios::scientific);
   std::cout.precision(6);

   std::cout << "Initializing flow field: " << std::endl;
   std::cout << "========================" << std::endl;
   for(i = 0; i < neqn+nauxvars; i ++){
     std::cout << "\tQinf[" << i << "] = " << this->Qinf[i] 
	       << "\t" << this->idata->GetDofName(i) << std::endl;
   }
   std::cout << std::endl;

   //set all the nodes interior and phantom
   for(i = 0; i < (nnode+gnode+nbnode); i++){
     for(j = 0; j < neqn; j++){
       this->space->q[i*(neqn+nauxvars) + j] = this->Qinf[j];
     }
     ComputeAuxiliaryVariables(&this->space->q[i*(neqn+nauxvars)]);
   }
}

template <class Type>
void HeatTransferEqnSet<Type>::ComputeAuxiliaryVariables(Type* Q)
{
  //do nothing, no aux variables
  return;
}

template <class Type>
Type HeatTransferEqnSet<Type>::MaxEigenvalue(Type* Q, Type* avec, Type vdotn, Type gamma, Type beta)
{
  return Q[0];
}


template <class Type>
void HeatTransferEqnSet<Type>::RoeFlux(Type* QL , Type* QR, Type* avec, Type vdotn, Type gamma, Type* flux, Type beta)
{
  //do nothing, we only have second order flux terms here
  flux[0] = 0.0;
}

template <class Type>
void HeatTransferEqnSet<Type>::ViscousFlux(Type* Q, Type* grad, Type* avec, Type mut, Type* flux)
{
  Type Tx = grad[0];
  Type Ty = grad[1];
  Type Tz = grad[2];
  
  //thermal conductivity
  Type k = 1.0;

  flux[0] = avec[3]*k*(Tx*avec[0] + Ty*avec[1] + Tz*avec[2]);
}

template <class Type>
void HeatTransferEqnSet<Type>::ViscousJacobian(Type* QL, Type* QR, Type* dx, Type s2,
					       Type* avec, Type mut, Type* aL, Type* aR)
{
  Int neqn = this->neqn;
  Type k = 1.0;
  Type dxnx = dx[0]*avec[0];
  Type dyny = dx[1]*avec[1];
  Type dznz = dx[2]*avec[2];

  Type c1 = k*avec[3]/s2;

  aL[0] = c1*(dxnx + dyny + dznz);

  for(Int i = 0; i < neqn*neqn; i++){
    aR[i] = aL[i];
  }
}

template <class Type>
Int HeatTransferEqnSet<Type>::GetGradientsLocation(std::vector<Int>& gradientLoc)
{
  gradientLoc.clear();
  //temperature
  gradientLoc.push_back(0);
  return gradientLoc.size();
}

template <class Type>
void HeatTransferEqnSet<Type>::ExtrapolateVariables(Type* Qho, const Type* Q, const Type* dQedge,
						    const Type* gradQ, const Type* dx, const Type* limiter)
{
  Qho[0] = Q[0] + this->ExtrapolateCorrection(dQedge[0], &gradQ[0*3], dx)*limiter[0];
}

//These are actually symmetry conditions, just named funny
template <class Type>
void HeatTransferEqnSet<Type>::GetInviscidWallBoundaryVariables(Type* QL, Type* QR, Type* Qinf, Type* avec, Type vdotn, Type beta)
{
  for(Int i = 0; i < this->neqn; i++){
    QR[i] = QL[i];
  }
  ComputeAuxiliaryVariables(QL);
  ComputeAuxiliaryVariables(QR);
}


template <class Type>
void HeatTransferEqnSet<Type>::GetIsothermalBoundaryVariables(Type* QL, Type* QR, Type Twall)
{
  for(Int i = 0; i < this->neqn; i++){
    QR[i] = QL[i] = Twall;
  }
  ComputeAuxiliaryVariables(QL);
  ComputeAuxiliaryVariables(QR);
}

template <class Type>
void HeatTransferEqnSet<Type>::GetHeatFluxBoundaryVariables(Type* QL, Type* QR, Type* normalQ, Type* normaldx, Type flux)
{
  //todo: this assumes that the inner point is completely normal
  //      it should be generalized for cases where we can't find a normal and only use the directional
  //      derivative part to set the boundary value

  //heat flux is defined as q = -k*(dT/dx) = -k * grad(T) - units are W/m^2, and others
  // q(+) is a flux into the solid, q(-) is out of the solid

  Type k = 1.0;
  //normaldx points from the wall to the infield point nearest to normal
    Type dist = sqrt(normaldx[0]*normaldx[0] + normaldx[1]*normaldx[1] + normaldx[2]*normaldx[2]);
  
  Type Twall = flux*dist/k + normalQ[0];

  //this is a hard set BC, set the surface and boundary condition both to the appropriate value to
  //recover the correct heat flux into the solid
  for(Int i = 0; i < this->neqn; i++){
    QR[i] = QL[i] = Twall;
  }

  ComputeAuxiliaryVariables(QL);
  ComputeAuxiliaryVariables(QR);
}
