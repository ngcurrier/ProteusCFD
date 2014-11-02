#include "macros.h"
#include "bc.h"
#include "eqnset.h"
#include "mesh.h"
#include "geometry.h"
#include "parallel.h"
#include "driver.h"
#include "crs.h"
#include "octree.h"
#include "solutionSpace.h"
#include "exceptions.h"

/*********************************************************/
//
// The code that follows is for the octree based method
// of solving for the wall distance.. keep separated
//
/********************************************************/

template <class Type>
void SyncParallelPoint(PObj<Type>* p, Type* xyzLocal, Type** xyzParallel, Int nptsLocal, Int* nptsParallel)
{
  Int i;
  Int sum;

  Int np = p->GetNp();
  Int rank = p->GetRank();

  Int countsRecv[np];
  Int offsetsRecv[np];

  MPI_Status status;
  MPI_Datatype mpit;
  MPI_Datatype mpitint;

  MPI_Request* sRequest = new MPI_Request[np];
  MPI_Request* rRequest = new MPI_Request[np];

  Type test = 0.0;
  mpit = MPI_GetType(test);
  mpitint = MPI_GetType(nptsLocal);

  for(i = 0; i < np; i++){
    countsRecv[i] = offsetsRecv[i] = 0;
  }

  //send the size of the block to every process
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Isend(&nptsLocal, 1, mpitint, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv the counts from other processes
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Irecv(&countsRecv[i], 1, mpitint, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }

  //sum up the nodes to expect
  sum = 0;
  for(i = 0; i < np; i++){
    sum += countsRecv[i];
  }
  offsetsRecv[0] = 0;
  for(i = 1; i < np; i++){
    offsetsRecv[i] = offsetsRecv[i-1] + countsRecv[i-1];
  }
  
  //allocate the memory to store the xyz coords
  *nptsParallel = sum;
  *xyzParallel = new Type[(*nptsParallel)*3];
  if(*xyzParallel == NULL){
    Abort << "WARNING: Out of memory in SyncParallelPoint()";
  }

    //send the size of the block to every process
  for(i = 0; i < np; i++){
    if(i != rank && nptsLocal != 0){
      MPI_Isend(xyzLocal, nptsLocal*3, mpit, i, 0, MPI_COMM_WORLD, &sRequest[i]);
    }
  }
  //recv the counts from other processes
  for(i = 0; i < np; i++){
    if(i != rank && countsRecv[i] != 0){
      MPI_Irecv(&((*xyzParallel)[offsetsRecv[i]*3]), countsRecv[i]*3, mpit, i, 0, MPI_COMM_WORLD, &rRequest[i]);
    }
  }
  //clear waits
  for(i = 0; i < np; i++){
    if(i != rank){
      MPI_Wait(&sRequest[i], &status);
      MPI_Wait(&rRequest[i], &status);
    }
  }


  delete [] sRequest;
  delete [] rRequest;

  return;
}

template <class Type>
void ComputeWallDistOct(Type* wallDist, SolutionSpace<Type>* space)
{
  Int i;

  BoundaryConditions<Real>* bc = space->bc;
  Mesh<Type>* m = space->m;
  PObj<Type>* p = space->p;

  //xyz coords of viscous nodes
  Type* xyzTotal = NULL;
  Type* xyzParallel = NULL;
  Type* xyzLocal = NULL;
  Int nptsTotal = 0;
  Int nptsLocal = 0;
  Int nptsParallel = 0;

  //octree parameters
  Int maxDepth = 12;
  Int minNodes = 15;

  //find all the viscous points in the mesh which are local
  //this would normally be handled within the driver module
  //but concern regarding the threading there kept me away
  for(Int eid = 0; eid < m->nbedge+m->ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int bcNum = bc->bc_map[factag];
    Int bcId;
    if(factag == 0){
      bcId = ParallelBoundary;
    }
    else{
      bcId = bc->bc_applied[bcNum];
    }  
    if(bcId == NoSlip){
      nptsLocal++;
    }
  }

  //allocate memory for viscous points
  xyzLocal = new Type[nptsLocal*3];
  if(xyzLocal == NULL){
    Abort << "WARNING: Allocation failure in ComputeWallDistOct() - xyzLocal";
  }

  Int count = 0;
  for(Int eid = 0; eid < m->nbedge+m->ngedge; eid++){
    Int factag = m->bedges[eid].factag;
    Int bcNum = bc->bc_map[factag];
    Int bcId;
    if(factag == 0){
      bcId = ParallelBoundary;
    }
    else{
      bcId = bc->bc_applied[bcNum];
    }  
    if(bcId == NoSlip){
      //the left node is the node on the surface
      Int left_cv = m->bedges[eid].n[0];
      memcpy(&xyzLocal[count*3], &m->xyz[left_cv*3], sizeof(Type)*3);
      count++;
    }
  }

  //update the viscous points from other processes
  SyncParallelPoint(p, xyzLocal, &xyzParallel, nptsLocal, &nptsParallel);

  //reallocate an array of the total length
  nptsTotal = nptsParallel+nptsLocal;
  xyzTotal = new Type[nptsTotal*3];
  if(xyzLocal == NULL){
    Abort << "WARNING: Allocation failure in ComputeWallDistOct() - xyzTotal";
  }

  //copy them over
  memcpy(xyzTotal, xyzLocal, sizeof(Type)*3*nptsLocal);
  delete [] xyzLocal;
  memcpy(&xyzTotal[nptsLocal*3], xyzParallel, sizeof(Type)*3*nptsParallel);
  delete [] xyzParallel;

  Octree<Type> octree(xyzTotal, nptsTotal, maxDepth, minNodes);
  octree.Print(std::cout);

  //search for the closest bounding box at the lowest level
  //then do an exhaustive search for all points in the sector
  for(i = 0; i < m->nnode+m->gnode; i++){
    Int nodeId;
    wallDist[i] = octree.FindDistance(&m->xyz[i*3 + 0], nodeId);
  }

  delete [] xyzTotal;

  return;
}

/*********************************************************/
//
// The code that follows is for the Laplacian based method
// of solving for the wall distance.. keep separated
//
/********************************************************/

//we solve linear advection equation with phi=0 at walls of interest
//and dphi/dn = 0 at other boundaries
//int_A{(del phi dot dn)} = -int_V{dv}
//d = sqrt{del phi (dot) del phi + 2phi} - abs{del phi}

template <class Type>
Type ComputeWallDist(Type* wallDist, SolutionSpace<Type>* space)
{
  Int i, j;
  Int maxiter = 100;
  Int iter;
  Int nsgs = 30;
  Type Cdes = 1.5;
  Real tol = 1.0e-12;
  Int gradType = 0; //use LSQ gradients

  Mesh<Type>* m = space->m;
  PObj<Type>* p = space->p;

  Type* phi = new Type[m->nnode+m->gnode+m->nbnode];
  Type* grads = new Type[3*(m->nnode+m->gnode)];
  Type norm;
  Type dphi;
  Type convnorm;

  CRS<Type> crs;
  crs.Init(m->nnode, m->gnode, 1, m->ipsp, m->psp, m->p);

  //setup gradient struct
  Gradient<Type> grad(1, 1, NULL, phi, space, gradType, grads, 0);
  //setup custom struct for flux/grad passing
  PointerStruct<Type> pstruct(crs.b, grads, phi);

  //init phi to one?  better guess?
  Type dx = m->extentsMax[0] - m->extentsMin[0];
  Type dy = m->extentsMax[1] - m->extentsMin[1];
  Type dz = m->extentsMax[2] - m->extentsMin[2];
  Type maxd = MAX(dx, dy);
  maxd = MAX(maxd, dz);
  for(i = 0; i < m->nnode; i++){
    Type ldx = pow(m->vol[i], 1.0/3.0);
    //phi[i] = maxd;
    //phi[i] = Cdes * ldx;
    phi[i] = maxd * ldx;
  }

  Kernel<Type> KWallDist(Kernel_WallDist);
  Kernel<Type> BKWallDist(BKernel_WallDist);
  Kernel<Type> BKWallBCs(BKernel_WallDistBC);
  Kernel<Type> WallDistJac(Kernel_WallDist_Jac);
  Kernel<Type> WallDistJacDiag(Kernel_WallDist_Diag);
  Kernel<Type> BWallDistJacDiag(BKernel_WallDist_Diag);

  crs.BlankSystem();
  //compute the jacobians... these don't change
  Driver(space, WallDistJac, 1, (void*)&crs);
  Driver(space, WallDistJacDiag, 1, (void*)&crs);
  Bdriver(space, BWallDistJacDiag, 1, (void*)&crs);
  //set initial guess to zero
  crs.BlankX();
  //take the inverse of the diagonal, once
  //not needed if using GMRES
  crs.A->PrepareSGS();
  //crs.A->PrintMatrix();
  
  iter = 0;
  do{
    //reset norms
    dphi = 0.0;
    norm = 0.0;

    //zero fluxes, blank system memory
    crs.BlankB();

    //update BC's
    BdriverNoScatter(space, BKWallBCs, 1, (void*)&pstruct);
    //compute gradients
    grad.Compute();
    //make sure gradients are zeroed on Neumann walls
    BdriverNoScatter(space, BKWallBCs, 1, (void*)&pstruct);
    //sync gradients via parallel call
    p->UpdateGeneralVectors(grads, 1*3);
   
    //scatter necessary fluxes
    Driver(space, KWallDist, 1, (void*)&pstruct);
    Bdriver(space, BKWallDist, 1, (void*)&pstruct);

    //add in source term
    for(j = 0; j < m->nnode; j++){
      crs.b[j] -= m->vol[j];
    }

    //solve
    convnorm = crs.SGS(nsgs, NULL, NULL, NULL, 1);
    //convnorm = crs.GMRES(1, 30, 2, NULL, NULL, NULL);

    for(j = 0; j < m->nnode; j++){
      phi[j] += crs.x[j];
    }
    p->UpdateGeneralVectors(phi, 1);

    norm = ParallelL2Norm(p, crs.b, m->nnode);
    dphi = ParallelL2Norm(p, crs.x, m->nnode);

    iter++;
    
    std::cout << iter << " ||res||:" << norm << " ||dphi||: " << dphi << " ||Solve||: " << norm << std::endl;

  }while(iter <= maxiter && real(norm) > tol);

  //compute gradient one last time
  grad.Compute();
  //sync gradients via parallel call
  p->UpdateGeneralVectors(grads, 1*3);
  
  std::cout << "Iters: " << iter << " ||dPhi||: " << dphi << " ||R||: " 
	    << norm << " ||D-SGS||: " << convnorm << std::endl;


  //now, compute walldistance
  for(i = 0; i < m->nnode; i++){
    wallDist[i] = abs(sqrt(abs(DotProduct(&grads[i*3], &grads[i*3]) + 2.0*phi[i])) - 
      Magnitude(&grads[i*3]));
  }

#if 0
  for(i = 0; i < m->nnode; i++){
    std::cout << m->xyz[i*3+0] << " " << m->xyz[i*3+1] << " " 
	      << m->xyz[i*3+2] << " " << phi[i] << " " << wallDist[i] << std::endl;
  }
#endif

  delete [] phi;
  delete [] grads;

  return norm;
}


template <class Type>
void Kernel_WallDist(KERNEL_ARGS)
{
  Int i;
  
  Mesh<Type>* m = space->m;

  PointerStruct<Type>* pstruct = (PointerStruct<Type>*) custom;
  Type* flux = pstruct->A;
  Type* grad = pstruct->B;
  Type* phi = pstruct->C;

  Type* gradL = &grad[left_cv*3];
  Type* gradR = &grad[right_cv*3];
  Type* phiL = &phi[left_cv];
  Type* phiR = &phi[right_cv];
  Type* avgGrad = (Type*)alloca(sizeof(Type)*3);
  Type dx[3];

  for(i = 0; i < 3; i++){
    avgGrad[i] = (gradL[i] + gradR[i]) / 2.0;
  }

  //use directional derivative method for Gauss point
  Type s2 = 0.0;
  Type qdots, dq;
  Type* xL = &m->xyz[left_cv*3];
  Type* xR = &m->xyz[right_cv*3];
  for(i = 0; i < 3; i++){
    dx[i] = xR[i] - xL[i];
    s2 += dx[i]*dx[i];
  }
  qdots = DotProduct(dx, avgGrad);
  dq = (phiR[0] - phiL[0] - qdots)/s2;
  for(i = 0; i < 3; i++){
    avgGrad[i] += dq*dx[i];
  }

  tempR[0] = (avgGrad[0]*avec[0] + avgGrad[1]*avec[1] + avgGrad[2]*avec[2])*avec[3];

  tempL[0] = -tempR[0];

  *size = 1;
  *ptrL = &flux[left_cv*1];
  *ptrR = &flux[right_cv*1];
  
  return;
}

template <class Type>
void BKernel_WallDist(B_KERNEL_ARGS)
{
  //TODO: should we include boundary edges in this computation?
  //      don't for now, we do include ghost(parallel) nodes though
  Int i;

  Mesh<Type>* m = space->m;

  PointerStruct<Type>* pstruct = (PointerStruct<Type>*) custom;
  Type* flux = pstruct->A;
  Type* grad = pstruct->B;
  Type* phi = pstruct->C;
  
  Type* gradL = &grad[left_cv*3];
  Type* gradR = &grad[right_cv*3];
  Type* phiL = &phi[left_cv];
  Type* phiR = &phi[right_cv];
  Type* avgGrad = (Type*)alloca(sizeof(Type)*3);
  Type* dx = (Type*)alloca(sizeof(Type)*3);
  Type* xL = &m->xyz[left_cv*3];

  if(m->IsGhostNode(right_cv)){
    Type* xR = &m->xyz[right_cv*3];
    for(i = 0; i < 3; i++){
      avgGrad[i] = (gradL[i] + gradR[i]) / 2.0;
    }
    //use directional derivative method for Gauss point
    Type s2 = 0.0;
    Type qdots, dq;
    for(i = 0; i < 3; i++){
      dx[i] = xR[i] - xL[i];
      s2 += dx[i]*dx[i];
    }
    qdots = DotProduct(dx, avgGrad);
    dq = (phiR[0] - phiL[0] - qdots)/s2;
    for(i = 0; i < 3; i++){
      avgGrad[i] += dq*dx[i];
    }
  }
  //we don't have grad information stored on BC's.. fake it
  else{
    for(i = 0; i < 3; i++){
      avgGrad[i] = gradL[i];
    }
  }


  tempL[0] = -(avgGrad[0]*avec[0] + avgGrad[1]*avec[1] + avgGrad[2]*avec[2])*avec[3];
  
  *size = 1;
  *ptrL = &flux[left_cv*1];
   
  return;
}


template <class Type>
PointerStruct<Type>::PointerStruct(Type* A_, Type* B_, Type* C_)
{
  A = A_;
  B = B_;
  C = C_;

  return;
}

template <class Type>
PointerStruct<Type>::~PointerStruct()
{
  return;
}


template <class Type>
void BKernel_WallDistBC(B_KERNEL_ARGS)
{
  //TODO: set only viscous walls to zero... set all others to zero derivative condition
  Int i;
  BoundaryConditions<Real>* bc = space->bc;
  Int bcNum = bc->bc_map[factag];
  Int bcId; 
  PointerStruct<Type>* pstruct = (PointerStruct<Type>*) custom;
  Type* flux = pstruct->A;
  Type* grad = pstruct->B;
  Type* phi = pstruct->C;

  if(factag == 0){
    bcId = ParallelBoundary;
  }
  else{
    bcId = bc->bc_applied[bcNum];
  }  
  if(bcId == NoSlip){
    phi[left_cv] = 0.0;
    phi[right_cv] = 0.0;
  }
  else if(bcId == ParallelBoundary){
    //do nothing.. updates take care of this
  }
  else{
    //make sure gradients are zero in normal direction
    Type dot = DotProduct(&grad[left_cv*3], avec);
    for(i = 0; i < 3; i++){
      grad[left_cv*3 + i] -= dot*avec[i];
    }
    phi[right_cv] = phi[left_cv];
  }
  *ptrL = NULL;
  *size = 0;
  
  return;
}

template <class Type>
void Kernel_WallDist_Jac(KERNEL_ARGS)
{
  Int i;

  Mesh<Type>* m = space->m;
  CRS<Type>* crs = (CRS<Type>*)custom;
  Type area = avec[3];

  Type ds2;
  Type* de = (Type*)alloca(sizeof(Type)*3);
  
  ds2 = 0.0;
  for(i = 0; i < 3; i++){
    de[i] = m->xyz[right_cv*3 + i] - m->xyz[left_cv*3 + i];
    ds2 += de[i]*de[i];
  }
  
  Type dx, dy, dz, d;
  dx = de[0]*avec[0];
  dy = de[1]*avec[1];
  dz = de[2]*avec[2];

  d = dx + dy + dz;
  
  //set data necessary for driver scatter
  *size = 1;
  *ptrL = crs->A->GetPointer(left_cv, right_cv);
  *ptrR = crs->A->GetPointer(right_cv, left_cv);

  tempL[0] = d*area/ds2;
  tempR[0] = -tempL[0];

  return;
}

template <class Type>
void Kernel_WallDist_Diag(KERNEL_ARGS)
{
  Mesh<Type>* m = space->m;
  CRS<Type>* crs = (CRS<Type>*) custom;

  //set data necessary for driver scatter
  *size = 1;
  *ptrL = crs->A->GetPointer(left_cv, left_cv);
  *ptrR = crs->A->GetPointer(right_cv, right_cv);

  Type* jacL = crs->A->GetPointer(left_cv, right_cv);
  Type* jacR = crs->A->GetPointer(right_cv, left_cv);

  tempL[0] = -jacL[0];
  tempR[0] = jacR[0];

  return;
}

template <class Type>
void BKernel_WallDist_Diag(B_KERNEL_ARGS)
{
  Int i;

  Mesh<Type>* m = space->m;
  CRS<Type>* crs = (CRS<Type>*)custom;
  Type area = avec[3];

  if(m->IsGhostNode(right_cv)){
    Type ds2;
    Type* de = (Type*)alloca(sizeof(Type)*3);
    ds2 = 0.0;
    for(i = 0; i < 3; i++){
      de[i] = m->xyz[right_cv*3 + i] - m->xyz[left_cv*3 + i];
      ds2 += de[i]*de[i];
    }
    
    Type dx, dy, dz, d;
    dx = de[0]*avec[0];
    dy = de[1]*avec[1];
    dz = de[2]*avec[2];
    
    d = dx + dy + dz;
  
    //set data necessary for driver scatter
    *size = 1;
    *ptrL = crs->A->GetPointer(left_cv, left_cv);
    
    tempL[0] = -d*area/ds2;
  }
  else{
    *size = 0;
    *ptrL = NULL;  
  }


  return;
}
