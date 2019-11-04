#include "bc_defines.h"
#include "bc.h"
#include "mesh.h"
#include "eqnset.h"
#include "driver.h"
#include "geometry.h"
#include "solutionSpace.h"

template <class Type>
Gradient<Type>::Gradient(Int nterms, Int stride, Int* list, Type* data,  
			 SolutionSpace<Type>* space, Int gradType, Type* grad, Bool weighted) :
  nterms(nterms), nnode(space->m->GetNumNodes() + space->m->GetNumParallelNodes()),  stride(stride), 
  type(gradType), weighted(weighted), space(space), data(data), grad(grad), list(NULL)
{
  if(data == NULL){
    std::cerr << "WARNING: In gradient constructor, data pointer is NULL!" << std::endl;
  }

  if(type == 1){
    std::cerr << "WARNING: Green Gauss gradients not yet tested as functional" << std::endl;
  }

  //assume we want all the values from zero up to nterms
  //unless something meaningful is passed in
  this->list = new Int[nterms];
  if(list == NULL){
    for(Int i = 0; i < nterms; i++){
      this->list[i] = i;
    }
  }
  else{
    for(Int i = 0; i < nterms; i++){
      this->list[i] = list[i];
    }
  }
  
  if(grad == NULL){
    //should be of size nterms*3*nnode
    Abort << "Gradient::Gradient() gradient memory is NULL external pointer";
  }
  if(data == NULL){
    Abort << "Gradient::Gradient() data memory is NULL external pointer";
  }
  if(space == NULL){
    Abort << "Gradient::Gradient() solution space memory is NULL external pointer";
  }
  
}

template <class Type>
Gradient<Type>::~Gradient()
{
  delete [] list;
}

template <class Type>
void Gradient<Type>::Compute()
{
  Int i; 
  Mesh<Type>* m = space->m;
  BoundaryConditions<Real>* bc = space->bc;

  //zero gradients
  for(i = 0; i < (nnode*nterms*3); i++){
    grad[i] = 0.0;
  }

  if(type == 0){
    //at this point node weighting coefficients should be
    //precomputed as they only change if the mesh moves
    Kernel<Type> Gradient(Kernel_LSQ_Gradient);
    Kernel<Type> BGradient(Bkernel_LSQ_Gradient);
    //call drivers to calculate gradients
    Driver(space, Gradient, nterms*3, this);
    Bdriver(space, BGradient, nterms*3, this);
  }
  else if(type == 1){
    Kernel<Type> Gradient(Kernel_Green_Gauss_Gradient);
    Kernel<Type> BGradient(Bkernel_Green_Gauss_Gradient);
    //call drivers to calculate gradients
    Driver(space, Gradient, nterms*3, this);
    Bdriver(space, BGradient, nterms*3, this);
    //divide gradients by the volume of the cell
    for(Int i = 0; i < nnode; ++i){
      for(Int j = 0; j < nterms*3; j++){
	grad[i*nterms*3 + j] /= m->vol[i];
      }
    }
  }
  //this ensures that the gradients are zero in all symmetry directions
  //we'll need a valid bc object to do this though, otherwise, badness
  if(bc != NULL){
    Kernel<Type> SymFix(Bkernel_Symmetry_Fix);
    BdriverNoScatter(space, SymFix, nterms*3, this);
  }

  //sync gradients via parallel call
  space->p->UpdateGeneralVectors(grad, nterms*3);

#if 0
  for(i = 0; i < nnode; i++){
    std::cout << "Gradient node: " << i << std::endl;
    for(Int j = 0; j < nterms; j++){
      std::cout << "Eqn: " << j << " " ;
      for(Int k = 0; k < 3; k++){
	std::cout << grad[i*nterms*3 + j*3 + k] << " " ;
      }
      std::cout << std::endl;
    }
  }
#endif
}

template <class Type>
void ComputeNodeLSQCoefficients(SolutionSpace<Type>* space)
{
  Mesh<Type>* m = space->m;
  Kernel<Type> LSQCoeff(Kernel_LSQ_Coefficients);
  Kernel<Type> BLSQCoeff(Bkernel_LSQ_Coefficients);
  Kernel<Type> LSQCoeffW(Kernel_LSQ_CoefficientsWeighted);
  Kernel<Type> BLSQCoeffW(Bkernel_LSQ_CoefficientsWeighted);

  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  for(Int i = 0; i < (nnode+gnode)*6; i++){
    m->s[i] = m->sw[i] = 0.0;
  }

  //call drivers to calculate coefficients
  Driver(space, LSQCoeff, 6, NULL);
  Bdriver(space, BLSQCoeff, 6, NULL);
  Driver(space, LSQCoeffW, 6, NULL);
  Bdriver(space, BLSQCoeffW, 6, NULL);
  
  //sync LSQ coefficients via parallel call
  space->p->UpdateGeneralVectors(space->m->s, 6);
  space->p->UpdateGeneralVectors(space->m->sw, 6);
}

template <class Type>
void ComputeLSQCoefficients(Type* s, Type* dxbar, Type* we)
{
  //expanded for clarity
  Type s11 = s[0];
  Type s12 = s[1];
  Type s13 = s[2];
  Type s22 = s[3];
  Type s23 = s[4];
  Type s33 = s[5];
  
  Type r11 = s11;
  Type r12 = s12;
  Type r13 = s13;
  Type r12_r11 = (real(r11) == 0.0) ? 0.0 : r12/r11;
  Type r22 = s22 - r12*r12_r11;
  Type r23 = s23 - r12_r11*r13;
  Type r13_r11 = (real(r11) == 0.0) ? 0.0 : r13/r11;
  Type r23_r22 = (real(r22) == 0.0) ? 0.0 : r23/r22;
  Type r33 = s33 - r13*r13_r11 - r23*r23_r22;
  Type dykdx = (dxbar[1] - (r12_r11)*dxbar[0]);

  //compute weight we_z
  we[2] = (real(r33) == 0.0) ? 0.0 : (dxbar[2] - r13_r11*dxbar[0] - r23_r22*dykdx)/r33;
  //compute weight we_y
  we[1] = (real(r22) == 0.0) ? 0.0 : (dykdx - r23*we[2])/r22;
  //compute weight we_x
  we[0] = (real(r11) == 0.0) ? 0.0 : (dxbar[0] - r12*we[1] - r13*we[2])/r11;
}


template <class Type>
void Kernel_Green_Gauss_Gradient(KERNEL_ARGS)
{
  Int i, j, k;
  
  Mesh<Type>* m = space->m;
  Gradient<Type>* g = (Gradient<Type>*) custom;

  Int nterms = g->GetNterms();
  Int nvars = g->GetStride();
  Int* list = g->list; 
  Type* QL = &g->data[left_cv*nvars];
  Type* QR = &g->data[right_cv*nvars];
  Type edgept[3];
  Type* faceavg = (Type*)alloca(sizeof(Type)*nterms);
  Type area = avec[3];

  //compute centroid of edge.. i.e. face centroid
  Centroid(&m->xyz[left_cv*3], &m->xyz[right_cv*3], edgept);

  //set data necessary for driver scatter
  *size = nterms*3;
  *ptrL = &g->grad[left_cv*nterms*3];
  *ptrR = &g->grad[right_cv*nterms*3];

  //compute face averaged q variables
  for(i = 0; i < nterms; i++){
    k = list[i];
    faceavg[i] = 0.5*(QL[k] +QR[k]);
  }
  
  for(i = 0; i < 3; i++){
    for(j = 0; j < nterms; j++){
      tempL[3*j + i] = faceavg[j]*avec[i]*area;
      tempR[3*j + i] = -faceavg[j]*avec[i]*area;
    }
  }
}

template <class Type>
void Bkernel_Green_Gauss_Gradient(B_KERNEL_ARGS)
{
  Int i, j, k;
  Mesh<Type>* m = space->m;
  Gradient<Type>* g = (Gradient<Type>*) custom;

  Int nterms = g->GetNterms();
  Int nvars = g->GetStride();
  Int* list = g->list;
  Type* QL = &g->data[left_cv*nvars];
  Type* QR = &g->data[right_cv*nvars];
  Type edgept[3];
  Type* faceavg = (Type*)alloca(sizeof(Type)*nterms);
  Type area = avec[3];

  //compute centroid of edge.. i.e. face centroid
  //TODO: fix this to use a point on the boundary itself.. this is wrong
  edgept[0] = m->xyz[left_cv*3 + 0];
  edgept[1] = m->xyz[left_cv*3 + 1];
  edgept[2] = m->xyz[left_cv*3 + 2];


  //set data necessary for driver scatter
  *size = nterms*3;
  *ptrL = &g->grad[left_cv*nterms*3];

  //compute face averaged q variables
  for(i = 0; i < nterms; i++){
    k = list[i];
    faceavg[i] = 0.5*(QL[k] + QR[k]);
  }
  
  for(i = 0; i < 3; i++){
    for(j = 0; j < nterms; j++){
      tempL[3*j + i] = faceavg[j]*avec[i]*area;
    }
  }
}

template <class Type>
void Kernel_LSQ_Gradient(KERNEL_ARGS)
{
  Int i, j, k;
  
  Mesh<Type>* m = space->m;
  Gradient<Type>* g = (Gradient<Type>*) custom;

  Int nterms = g->GetNterms();
  Int nvars = g->GetStride();
  Int* list = g->list;
  Type dx[3];
  Type* xL = &m->xyz[left_cv*3];
  Type* xR = &m->xyz[right_cv*3];
  
  Subtract(xL, xR, dx);
  
  Type* qL = &g->data[left_cv*nvars];
  Type* qR = &g->data[right_cv*nvars];

  Type* sL = &m->s[left_cv*6];
  Type* sR = &m->s[right_cv*6];
  if(g->IsWeighted()){
    sL = &m->sw[left_cv*6];
    sR = &m->sw[right_cv*6];
  }

  Type *gradL = &g->grad[left_cv*nterms*3];
  Type *gradR = &g->grad[right_cv*nterms*3];

  Type weL[3];
  Type weR[3];
  Type dq;
  Type weight = 1.0;

  if(g->IsWeighted()){
    //use inverse distance weighting
    Type dx2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    weight = 1.0/sqrt(dx2);
    dx[0] *= weight;
    dx[1] *= weight;
    dx[2] *= weight;
  }
  //compute edge coefficients for left_cv
  ComputeLSQCoefficients(sL, dx, weL);
  
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];

  //compute edge coefficients for right_cv
  ComputeLSQCoefficients(sR, dx, weR);

  //set data necessary for driver scatter
  *size = nterms*3;
  *ptrL = gradL;
  *ptrR = gradR;

  for(i = 0; i < nterms; i++){
    k = list[i];
    dq = weight*(qR[k] - qL[k]);
    for(j = 0; j < 3; j++){
      tempL[3*i + j] = -weL[j]*dq;
      tempR[3*i + j] = +weR[j]*dq;
    }
  }
}

template <class Type>
void Bkernel_LSQ_Gradient(B_KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Int i, j, k;
    Gradient<Type>* g = (Gradient<Type>*) custom;
    Int nterms = g->GetNterms();
    Int nvars = g->GetStride();
    Int* list = g->list;
    Type dx[3];
    Type* xL = &m->xyz[left_cv*3];
    Type* xR = &m->xyz[right_cv*3];
    
    Subtract(xL, xR, dx);
    
    Type* qL = &g->data[left_cv*nvars];
    Type* qR = &g->data[right_cv*nvars];
    
    Type* sL = &m->s[left_cv*6];
    if(g->IsWeighted()){
      sL = &m->sw[left_cv*6];
    }

    Type *gradL = &g->grad[left_cv*nterms*3];
    
    Type weL[3];
    Type dq;
    Type weight = 1.0;
    
    if(g->IsWeighted()){
      //use inverse distance weighting
      Type dx2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
      weight = 1.0/sqrt(dx2);
      dx[0] *= weight;
      dx[1] *= weight;
      dx[2] *= weight;
    }
    //compute edge coefficients for left_cv
    ComputeLSQCoefficients(sL, dx, weL);
    
    //set data necessary for driver scatter
    *size = nterms*3;
    *ptrL = gradL;
    *ptrR = NULL;
    
    for(i = 0; i < nterms; i++){
      k = list[i];
      dq = weight*(qR[k] - qL[k]);
      for(j = 0; j < 3; j++){
	tempL[3*i + j] = -weL[j]*dq;
      }
    }
  }
  else{
    *size = 0;
    *ptrL = NULL;
    *ptrR = NULL;
  }
}

template <class Type>
void Kernel_LSQ_Coefficients(KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;

  Type dx[3];
  Type* x1 = &m->xyz[left_cv*3];
  Type* x2 = &m->xyz[right_cv*3];
  
  Subtract(x1, x2, dx);
  
  Type* s1 = &m->s[left_cv*6];
  Type* s2 = &m->s[right_cv*6];

  Type dx2 = dx[0]*dx[0];
  Type dy2 = dx[1]*dx[1];
  Type dz2 = dx[2]*dx[2];
  
  //add sum from this edge, s1
  tempL[0] = dx2;             //s11
  tempL[1] = (dx[0]*dx[1]);   //s12
  tempL[2] = (dx[0]*dx[2]);   //s13
  tempL[3] = dy2;             //s22
  tempL[4] = (dx[1]*dx[2]);   //s23
  tempL[5] = dz2;             //s33

  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];

  //squared values are always positive, no need to be cautious of signs
  tempR[0] = dx2;             //s11
  tempR[1] = (dx[0]*dx[1]);   //s12
  tempR[2] = (dx[0]*dx[2]);   //s13
  tempR[3] = dy2;             //s22
  tempR[4] = (dx[1]*dx[2]);   //s23
  tempR[5] = dz2;             //s33

  *size = 6;  
  *ptrL = s1;
  *ptrR = s2;
}

template <class Type>
void Bkernel_LSQ_Coefficients(B_KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Type dx[3];
    Type* x1 = &m->xyz[left_cv*3];
    Type* x2 = &m->xyz[right_cv*3];
    
    Subtract(x1, x2, dx);
    
    Type* s1 = &m->s[left_cv*6];
    
    Type dx2 = dx[0]*dx[0];
    Type dy2 = dx[1]*dx[1];
    Type dz2 = dx[2]*dx[2];
    
    //add sum from this edge, s1
    tempL[0] = dx2;                   //s11
    tempL[1] = (dx[0]*dx[1]);         //s12
    tempL[2] = (dx[0]*dx[2]);         //s13
    tempL[3] = dy2;                   //s22
    tempL[4] = (dx[1]*dx[2]);         //s23
    tempL[5] = dz2;                   //s33
    
    //these are going to be weighted later so we're done here
    
    *size = 6;  
    *ptrL = s1;
  }
  else{
    *size = 0;
    *ptrL = NULL;
  }
}

template <class Type>
void Kernel_LSQ_CoefficientsWeighted(KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;

  Type dx[3];
  Type* x1 = &m->xyz[left_cv*3];
  Type* x2 = &m->xyz[right_cv*3];
  
  Subtract(x1, x2, dx);
  
  Type* s1 = &m->sw[left_cv*6];
  Type* s2 = &m->sw[right_cv*6];

  Type ds2 = 1.0/(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

  Type dx2 = dx[0]*dx[0]*ds2;
  Type dy2 = dx[1]*dx[1]*ds2;
  Type dz2 = dx[2]*dx[2]*ds2;
  
  //add sum from this edge, s1
  tempL[0] = dx2;               //s11
  tempL[1] = dx[0]*dx[1]*ds2;   //s12
  tempL[2] = dx[0]*dx[2]*ds2;   //s13
  tempL[3] = dy2;               //s22
  tempL[4] = dx[1]*dx[2]*ds2;   //s23
  tempL[5] = dz2;               //s33

  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];

  //squared values are always positive, no need to be cautious of signs
  tempR[0] = dx2;                 //s11
  tempR[1] = (dx[0]*dx[1])*ds2;   //s12
  tempR[2] = (dx[0]*dx[2])*ds2;   //s13
  tempR[3] = dy2;                 //s22
  tempR[4] = (dx[1]*dx[2])*ds2;   //s23
  tempR[5] = dz2;                 //s33

  *size = 6;  
  *ptrL = s1;
  *ptrR = s2;
}

template <class Type>
void Bkernel_LSQ_CoefficientsWeighted(B_KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;

  if(m->IsGhostNode(right_cv)){
    Type dx[3];
    Type* x1 = &m->xyz[left_cv*3];
    Type* x2 = &m->xyz[right_cv*3];
    
    Subtract(x1, x2, dx);
    
    Type* s1 = &m->sw[left_cv*6];
    
    Type ds2 = 1.0/(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    
    Type dx2 = dx[0]*dx[0]*ds2;
    Type dy2 = dx[1]*dx[1]*ds2;
    Type dz2 = dx[2]*dx[2]*ds2;
    
    //add sum from this edge, s1
    tempL[0] = dx2;
    tempL[1] = dx[0]*dx[1]*ds2;
    tempL[2] = dx[0]*dx[2]*ds2;
    tempL[3] = dy2;
    tempL[4] = dx[1]*dx[2]*ds2;
    tempL[5] = dz2;
    
    //these are going to be weighted later so we're done here
    
    *size = 6;  
    *ptrL = s1;
  }
  else{
    *size = 0;
    *ptrL = NULL;
  }
}

template <class Type>
void Bkernel_Symmetry_Fix(B_KERNEL_ARGS)
{
  Int i, j;
  Gradient<Type>* g = (Gradient<Type>*) custom;
  BoundaryConditions<Real>* bc = space->bc;
  Mesh<Type>* m = space->m;
  Int bcType = space->bc->GetBCType(factag); 
  Int nterms = g->GetNterms();
  Type dot;
  Type* gradL = &g->grad[left_cv*nterms*3];

  if(bcType == Proteus_Symmetry){
    //ensure gradients are zero in the normal direction
    for(i = 0; i < nterms; i++){
      dot = DotProduct(&gradL[i*3], avec);
      for(j = 0; j < 3; j++){
	gradL[i*3 + j] -= dot*avec[j];
      }
    }
  }
}

template <class Type>
Type DirectionalDerivativeGrad(Type* dx, Type s2, Type* avec, Type* qL, Type* qR, Type* grad, 
			       std::vector<Int> gradLoc)
{
  Int i, j;
  
  Type dq;
  for(j = 0; j < gradLoc.size(); j++){
    Int loc = gradLoc[j];
    dq = (qR[loc] - qL[loc])/s2;
    for(i = 0; i < 3; i++){
      grad[j*3 + i] += dq*dx[i];
    }
  }

  //this is the approximate derivative of the gradient
  return (-DotProduct(dx, avec)/s2);
}
