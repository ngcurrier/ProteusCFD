#include "macros.h"
#include "bc.h"
#include "matrix.h"
#include "geometry.h"
#include "mesh.h"
#include "mem_util.h"
#include "parallel.h"
#include "parametric.h"
#include "walldist.h"
#include "pythonInterface.h"
#include <cmath>

//these defines are for the sinusoidally oscillating airfoil case
#define ALPHA_VARY 2.41 //degrees
#define ALPHA_MEAN 2.89 //degrees
#define OSC_FREQ 50.32  //hertz

template <class Type>
void MoveMesh(SolutionSpace<Type>* space, void* custom){
  Int i;
  Mesh<Type>* m = space->m;
  Int nnode = m->GetNumNodes();
  Param<Type>* param = space->param;
  BoundaryConditions<Real>* bc = space->bc;
  Type* dx = new Type[nnode*3];

  //default to not requiring mesh smoothing for dx values
  Int reqSmoothing = false;

  //call routine to get point displacements on boundaries
  MemBlank(dx, nnode*3);
  UpdateMeshDisplacements(space, dx, &reqSmoothing, custom);

  if(reqSmoothing){
    //use linear elastic to move the mesh
    Int smoothingPasses = 100;
    MoveMeshLinearElastic(m, bc, dx, smoothingPasses);
  }
  else{
    //dx values are explicit, add them up
    for(i = 0; i < nnode*3; i++){
      m->xyz[i] += dx[i];
    }
    m->p->UpdateXYZ(m->xyz);
  }

  //update mesh metrics
  m->CalcMetrics();

  //update mesh velocities using total motion over all the Newton steps
  UpdateMeshVelocities(m, m->xyz, m->xyzold, m->xyzoldm1, param->dt, space->iter, param->torder);

  //update gradient weighting coefficients
  ComputeNodeLSQCoefficients(space);

  delete [] dx;

  return;
}

template <class Type>
void UpdateMeshDisplacements(SolutionSpace<Type>* space, Type* dx, Int* reqSmoothing, void* custom)
{
  Int i;
  
  Mesh<Type>* m = space->m;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = space->param;

  Int nnode = m->GetNumNodes();

  //sinusoidally pitching airfoil (rigid body)
  //this only works if we are rotating about the origin
  if(param->movement == 1){
    Type pi = acos(-1.0);
    //number of rads to move both up and down
    Type alphaVaryRad = ALPHA_VARY*pi/180.0;
    //the mean angle of attack
    Type alphaMeanRad = ALPHA_MEAN*pi/180.0;
    //get current position of airfoil
    Type dimOmega = OSC_FREQ*2.0*pi;      //(rad/s)
    //non-dimensional circular frequency
    Type nondimOmega = dimOmega*param->ref_time;
    //compute omega*t
    Type radPos = nondimOmega*(Type)space->iter*param->dt;
    //since sine varies from 0 -> 1 -> 0 -> -1 over 2*pi
    //we can simply multiply by our alpha variation to get
    //a smooth pitching action
    Type alphaRad = sin(radPos)*alphaVaryRad;

    //we use the negative b/c of the way a rotation matrix works
    //by default
    Type sn = sin(-alphaRad);
    Type cs = cos(-alphaRad);

    Type nxyz[3];
    //compute the displacements for the interior nodes
    for(i = 0; i < nnode; i++){
      //use the mesh at time t=0 for the offset
      Type* pt_base = &m->xyz_base[i*3];
      Type* pt_now = &m->xyz[i*3];

      //apply the rotation matrix
      // | cos(theta) -sin(theta) |
      // | sin(theta)  cos(theta) |
      // this matrix by default rotates CCW for positive theta, we use the negative
      // of this below b/c it is the opposite

      nxyz[0] = cs*pt_base[0] - sn*pt_base[1];
      nxyz[1] = sn*pt_base[0] + cs*pt_base[1];

      dx[i*3 + 0] = nxyz[0] - pt_now[0];
      dx[i*3 + 1] = nxyz[1] - pt_now[1];
      //we are rotating in the x-y plane only
      dx[i*3 + 2] = 0.0;
    }
#if 0
    //for this particular case we also have to update the lift
    //direction that was read from the param file since we are rotating
    //the whole mesh...here we are actually measuring CN coeff not lift
    Type alphaLift = alphaRad;
    param->liftdir[0] = real(sin(alphaLift));
    param->liftdir[1] = real(cos(alphaLift));
    param->liftdir[2] = 0.0;
#endif
    //do not require that a mesh movement routine be called, all
    //displacements are fully determined here
    *reqSmoothing = false;
  }
  //sinusoidally pitching airfoil (moving grid)
  //this only works if we are rotating about the origin
  else if(param->movement == 2){
    Type pi = acos(-1.0);
    //number of rads to move both up and down
    Type alphaVaryRad = ALPHA_VARY*pi/180.0;
    //the mean angle of attack
    Type alphaMeanRad = ALPHA_MEAN*pi/180.0;
    //get current position of airfoil
    Type dimOmega = OSC_FREQ*2.0*pi;      //(rad/s)
    //non-dimensional circular frequency
    Type nondimOmega = dimOmega*param->ref_time;
    //compute omega*t
    Type radPos = nondimOmega*(Type)space->iter*param->dt;
    //since sine varies from 0 -> 1 -> 0 -> -1 over 2*pi
    //we can simply multiply by our alpha variation to get
    //a smooth pitching action
    Type alphaRad = sin(radPos)*alphaVaryRad;

    //we use the negative b/c of the way a rotation matrix works
    //by default
    Type sn = sin(-alphaRad);
    Type cs = cos(-alphaRad);

    Type nxyz[3];

    //retrieve that nodes that are on the design surface indicated in the
    //boundary conditions file
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, bc, &nodelist);
    //compute the displacements for the interior nodes
    for(i = 0; i < nodecount; i++){
      Int node = nodelist[i];
      //use the mesh at time t=0 for the offset
      Type* pt_base = &m->xyz_base[node*3];
      Type* pt_now = &m->xyz[node*3];
      //apply the rotation matrix
      // | cos(theta) -sin(theta) |
      // | sin(theta)  cos(theta) |
      // this matrix by default rotates CCW for positive theta, we use the negative
      // of this below b/c it is the opposite

      nxyz[0] = cs*pt_base[0] - sn*pt_base[1];
      nxyz[1] = sn*pt_base[0] + cs*pt_base[1];

      dx[node*3 + 0] = nxyz[0] - pt_now[0];
      dx[node*3 + 1] = nxyz[1] - pt_now[1];
      //we are rotating in the x-y plane only
      dx[node*3 + 2] = 0.0;
    }
#if 0
    //for this particular case we also have to update the lift
    //direction that was read from the param file since we are rotating
    //the airfoil...here we are actually measuring CN coeff not lift
    Type alphaLift = alphaRad;
    param->liftdir[0] = real(sin(alphaLift));
    param->liftdir[1] = real(cos(alphaLift));
    param->liftdir[2] = 0.0;
#endif

    //clear memory for node list
    delete [] nodelist;

    //require that a mesh movement routine be called, only
    //displacements given are on the surface
    *reqSmoothing = true;
  }
  //this is the movement routine for specified pitching and plunge, via
  //aeroelastic typical section routine
  else if(param->movement == 3 || param->movement == 6){
    ScalarField<Type>& pitch = space->GetScalarField("Pitch");
    ScalarField<Type>& plunge = space->GetScalarField("Plunge");
    //these are dimensional values, we need to nondimensionalize them for the
    //CFD side of things
    Type h = plunge.GetField()/param->ref_length;
    Type alpha = pitch.GetField();

    //we use the negative b/c of the way a rotation matrix works
    //by default
    Type sn = sin(-alpha);
    Type cs = cos(-alpha);

    Type nxyz[3];
    //retrieve that nodes that are on the design surface indicated in the
    //boundary conditions file
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, bc, &nodelist);

    //compute the displacements for the interior nodes
    for(i = 0; i < nodecount; i++){
      Int node = nodelist[i];
      //use the mesh at time t=0 for the offset
      Type* pt_base = &m->xyz_base[node*3];
      Type* pt_now = &m->xyz[node*3];

      //apply the rotation matrix
      // | cos(theta) -sin(theta) |
      // | sin(theta)  cos(theta) |
      // this matrix by default rotates CCW for positive theta, we use the negative
      // of this below b/c it is the opposite

      nxyz[0] = cs*pt_base[0] - sn*pt_base[1];
      //plunge assumed in the y-direction, positive down
      nxyz[1] = sn*pt_base[0] + cs*pt_base[1] - h;
      nxyz[2] = 0.0;

      dx[node*3 + 0] = nxyz[0] - pt_now[0];
      dx[node*3 + 1] = nxyz[1] - pt_now[1];
      //we are rotating in the x-y plane only
      dx[node*3 + 2] = 0.0;
    }
    //NOTE: for this particular case, we do not update the lift direction
    //since we assume for typical section airfoil analysis that the plunge
    //is calculated via Cl and the wing root is fixed... thus, any twisting
    //should not affect the force direction the structure cares about

    //clear memory for node list
    delete [] nodelist;

    //assume the elastic axis is centered at (0,0,0)
    //since the airfoil is plunging we have to move the elastic axis
    //definition with the body, use the first body defined in the BC file
    space->forces->bodies[1].momentPt[0] = 0.0;
    space->forces->bodies[1].momentPt[1] = -h;
    space->forces->bodies[1].momentPt[2] = 0.0;

    //require that a mesh movement routine be called, only
    //displacements given are on the surface
    *reqSmoothing = true;
  }
  else if(param->movement == 4){
    //this is a testing routine, move at speed 0.6 mach in x direction
    for(i = 0; i < nnode; i++){
      dx[i*3 + 0] = 0.6*param->dt;
      dx[i*3 + 1] = 0.0;
      dx[i*3 + 2] = 0.0;
    }
    *reqSmoothing = false;
  }
  else if(param->movement == 5){
    //this is a testing routine, vibrate vertically 0.1
    //retrieve that nodes that are on the design surface indicated in the
    //boundary conditions file
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, bc, &nodelist);
    //compute the displacements for the interior nodes
    for(i = 0; i < nodecount; i++){
     Int node = nodelist[i];
     //for(Int node = 0; node < nnode; node++){
     dx[node*3 + 0] = 0.0;
     dx[node*3 + 1] = 0.1*sin(1.5*(space->iter)*param->dt);
     dx[node*3 + 2] = 0.0;
    }
    *reqSmoothing = true;
  }
  else if(param->movement == 7){
    double* dxyz = (double*)custom;
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, bc, &nodelist);
    //compute the displacements for the interior nodes
    for(i = 0; i < nodecount; i++){
     Int node = nodelist[i];
     dx[node*3 + 0] = dxyz[i*3 + 0];
     dx[node*3 + 1] = dxyz[i*3 + 1];
     dx[node*3 + 2] = dxyz[i*3 + 2];
    }
    *reqSmoothing = true;
  }
  else if(param->movement == 8){ //python boundary movement calls
    double* dxyz = (double*)custom;
    Int* nodelist = NULL;
    Int nodecount = GetNodesOnMovingBC(m, bc, &nodelist);
#ifdef _HAS_PYTHON
    if (MPI_GetType(m->xyz_base[0]) == MPI_COMPLEX){
      Abort << "PythonInterface::GetBoundaryMovement() cannot handle complex variable types";
    }
    PythonWrapper pywrap("./", "getBoundaryMovement", "getBoundaryMovement");
    pywrap.GetBoundaryMovement(real(Type(space->iter)*param->dt), nodelist, nodecount, m->GetNumNodes(),
			       (double*)m->xyz_base, dxyz);
#else
    Abort << "UpdateMeshDisplacements() - python not built with solver";
#endif
    
    *reqSmoothing = true;
  }
}

template <class Type, class Type2>
void UpdateMeshVelocities(Mesh<Type>* m, Type* xyz, Type* xyzold, Type* xyzoldm1, 
			  Type2 timestep, Int iter, Int torder)
{
  Int i;
  Type cnm1, cn, cnp1;
  Type K0, K1, K2;

  if(iter > 2 && torder == 2){
    //coefficients for BDF2 are phi_n+1 = 1.5, phi_n = -2.0, phi_n-1 = 0.5 
    cnp1 = 1.5;
    cnm1 = -0.5;
    cn = -2.0;
  }
  else{
    cnp1 = 1.0;
    cn = -1.0;
    cnm1 = 0.0;
  }

  for(i = 0; i < m->GetNumNodes(); i++){
    //this is the BDF2 grid speed
    K0 = (cnp1*xyz[i*3 + 0] + cn*xyzold[i*3 + 0] - cnm1*xyzoldm1[i*3 + 0]);
    K1 = (cnp1*xyz[i*3 + 1] + cn*xyzold[i*3 + 1] - cnm1*xyzoldm1[i*3 + 1]);
    K2 = (cnp1*xyz[i*3 + 2] + cn*xyzold[i*3 + 2] - cnm1*xyzoldm1[i*3 + 2]);
    m->nv[i*3 + 0] = K0/timestep;
    m->nv[i*3 + 1] = K1/timestep;
    m->nv[i*3 + 2] = K2/timestep;
  }
  
  m->p->UpdateGeneralVectors(m->nv, 3);

  return;
}

template <class Type>
void SmoothMeshLaplacian(Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			 const Real tol, const Int max_iter, 
			 const bool no_movement_dir[3], const Int weighting)
{
  Int i,j,k,node,onode;
  Int indx1, indx2, iter;
  Type w;
  Type* dx = new Type[3];
  Type norm;
  Int normn;

  PObj<Type>* p = m->p;

  Int* symNodes;
  Int nSymNodes = GetNodesOnSymmetryPlanes(m, bc, &symNodes);

  iter = 0;

  do{
    //loop over interior nodes
    //and compute weights
    normn = 0;
    norm = 0.0;
    for(i = 0; i < m->nvnodes+nSymNodes; i++){
      if(i < m->nvnodes){
	node = m->vnodes[i];
      }
      else{
	node = symNodes[i - m->nvnodes];
      }
      indx1 = m->ipsp[node];
      indx2 = m->ipsp[node+1];
      dx[0] = dx[1] = dx[2] = 0.0;
      //simple weighting
      //scale dependent weighting
      for(j = indx1; j < indx2; j++){
	onode = m->psp[j];
	if(weighting == 0){
	  w = 1.0 / CAbs((Type)(indx2-indx1));
	}
	else if(weighting == 1){
	  w = 1.0 / Distance(&m->xyz[node*3 + 0], &m->xyz[onode*3 + 0]);
	}
	else{
	  std::cerr << "Weighting method " << weighting << " unrecognized" << std::endl;
	  return;
	}
	//update movement vector
	dx[0] += w * (m->xyz[onode*3 + 0] - m->xyz[node*3 + 0]);
	dx[1] += w * (m->xyz[onode*3 + 1] - m->xyz[node*3 + 1]);
	dx[2] += w * (m->xyz[onode*3 + 2] - m->xyz[node*3 + 2]);

      }
#if 0
      Type max_step = 1.0e-1;
      //limit maximum step
      for(k = 0; k < 3; k++){
	dx[k] = MIN(dx[k], max_step);
      }
#endif
      //update point and movement norm
      //move mesh points
      for(k = 0; k < 3; k++){
	if(!no_movement_dir[k]){
	  m->xyz[node*3 + k] += dx[k];
	  norm += (dx[k]*dx[k]);
	  normn++;
	}
      }
    }
    p->UpdateXYZ(m);
    iter++;
    norm = sqrt(norm) / (Type)normn;
    norm = ParallelL2Norm(p, &norm, 1);
    std::cout << "Laplace smoothing |dx|: " << norm << std::endl;
    
  }while(real(norm) > tol && iter < max_iter);
  std::cout << "Laplace smoothing converged after " << iter << " iterations" << std::endl;
  if(iter >= max_iter){
    std::cout << "Reached maximum number of iterations" << std::endl;
  }

  delete [] dx;
  delete [] symNodes;

  return;
}

template <class Type>
Int BumpNodesOnLine(const Type dist[3], const Int ptId, const Int nnodes, 
		    Type* xyz, const Int* ipsp,const Int* psp, 
		    const Int direction)
{
  Int i, j, indx1, indx2;
  Int node2;
  Type x[2];
  //search tolerance for line
  Real tol = 1.0e-8;
  Int count = 0;
  //TODO: is there a better starting point here?
  Int mem_size = 50;
  Int* list = new Int[mem_size];
  //set passed in node as first in the list
  list[0] = ptId;
  Bool* searched = new Bool[nnodes];

  for(i = 0; i < nnodes; i++){
    searched[i] = 0;
  }
  searched[ptId] = 1;

  Int dir1, dir2;
  if(direction == 0){
    dir1 = 1; 
    dir2 = 2;
    x[dir1] = xyz[ptId*3 + dir1];
    x[dir2] = xyz[ptId*3 + dir2];
  }
  else if(direction == 1){
    dir1 = 0;
    dir2 = 2;
    x[dir1] = xyz[ptId*3 + dir1];
    x[dir2] = xyz[ptId*3 + dir2];
  }
  else if(direction == 2){
    dir1 = 0; 
    dir2 = 1;
    x[dir1] = xyz[ptId*3 + dir1];
    x[dir2] = xyz[ptId*3 + dir2];
  }
  else{
    std::cerr << "Invalid direction! Pick 0 (x), 1(y), or 2(z)";
    return (-1);
  }

  i = 0;
  do{
    indx1 = ipsp[list[i]];
    indx2 = ipsp[list[i]+1];
    for(j = indx1; j < indx2; j++){
      node2 = psp[j];
      //if node is within tolerance of line in correct direction
      //increase counter
      if(real(Abs(xyz[node2*3 + dir1] - x[dir1]) < tol) && real(Abs(xyz[node2*3 + dir2] - x[dir2]) < tol) && !searched[node2]){
	//if the list runs out of memory, make it bigger
	searched[node2] = 1;
	if((count+1) >= mem_size){
	  MemResize(&list, mem_size, mem_size+50);
	  mem_size += 50;
	} 
	count++;
	list[count] = node2;
      }
    }
    i++;
  }while(i <= count);

  std::cout << "NODES BUMPED: " << std::endl;
  //bump the nodes
  for(j = 0; j <= count; j++){
    std::cout << list[j] << " ";
    node2 = list[j];
    xyz[node2*3 + 0] += dist[0];
    xyz[node2*3 + 1] += dist[1];
    xyz[node2*3 + 2] += dist[2];
  }
  std::cout << std::endl;

  delete [] list;
  delete [] searched;

  return (count+1);
}

template <class Type>
void MoveMeshLinearElastic(Mesh<Type>* m, BoundaryConditions<Real>* bc, Type* dx, const Int iter)
{
  //code follows AIAA 2005-0923 (Karman, et. al.)

  Int i, j, eid;
  Type deltaDq = 0.0;

  Type* ptr1;
  Type* ptr2;
  Type* ptrd1;
  Type* ptrd2;
  Type* avec;
  Type area;
  Type* aspectRatio;
  Type E; //young's modulus
  Type nu = 0.2; //poisson's ratio
  Type nup1 = 1.0 + nu;
  Type nup2 = 1.0 - 2.0*nu;
  Type fact1 = (1.0-nu)/(nup1*nup2);
  Type fact2 = nu/(nup1*nup2);
  Type ds2, dot;
  Type de[3];
  Int dim2 = 9;
  Int n1, n2;

  //alphas
  Type a11, a12, a13, a21, a22, a23, a31, a32, a33;
  //thetas
  Type t11, t12, t13, t21, t22, t23, t31, t32, t33;

  
  CRS<Type> crs;
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int nedge = m->GetNumEdges();
  Int ngedge = m->GetNumParallelEdges();
  Int nbedge = m->GetNumBoundaryEdges();
  crs.Init(nnode, gnode, 3, m->ipsp, m->psp, m->p);

  aspectRatio = new Type[nnode];
  for(i = 0; i < nnode; i++){
    aspectRatio[i] = 1.0/m->vol[i];
  }  

  //blank system memory
  crs.BlankSystem();

  //set the initial guess for the system
  memcpy(crs.x, dx, (nnode*3)*sizeof(Type));

  for(eid = 0; eid < nedge; eid++){
    avec = m->edges[eid].a;
    area = avec[3];
    n1 = m->edges[eid].n[0];
    n2 = m->edges[eid].n[1];
    E = (aspectRatio[n1] + aspectRatio[n2])/2.0;

    //This section based on directional derivatives as given
    //Hyam's disseration p. 21
    ds2 = 0.0;
    dot = 0.0;
    for(j = 0; j < 3; j++){
      de[j] = m->xyz[n2*3 + j] - m->xyz[n1*3 + j];
      ds2 += de[j]*de[j];
    }
    for(j = 0; j < 3; j++){
      dot += de[j]*avec[j];
    }
    dot *= area;

    a11 = a22 = a33 = E*fact1;
    a12 = a13 = a21 = a23 = a31 = a32 = 0.5*E/nup1;
    t11 = t22 = t33 = E*fact2;
    t12 = t13 = t21 = t23 = t31 = t32 = 0.5*E/nup1;
    
    //off diagonal entry for row i1, column j1 -- dfQ(node1)/dQ(node2)
    //normal vector points from i1 to i2, differentiate edge flux w.r.t. to QR
    ptr1 = crs.A->GetPointer(n1, n2);
    ptr2 = crs.A->GetPointer(n2, n1);
    ptr1[0] = (a11+a12+a13)*dot/ds2;
    ptr1[3] = (t11*de[1]*avec[0] + t12*de[0]*avec[1])*area/ds2;
    ptr1[6] = (t11*de[2]*avec[0] + t13*de[0]*avec[2])*area/ds2;
    
    ptr1[1] = (t21*de[1]*avec[0] + t22*de[0]*avec[1])*area/ds2;
    ptr1[4] = (a21+a22+a23)*dot/ds2; 
    ptr1[7] = (t22*de[2]*avec[1] + t23*de[1]*avec[2])*area/ds2;
    
    ptr1[2] = (t31*de[2]*avec[0] + t33*de[0]*avec[2])*area/ds2;
    ptr1[5] = (t32*de[2]*avec[1] + t33*de[1]*avec[2])*area/ds2;
    ptr1[8] = (a31+a32+a33)*dot/ds2;	       

    //copy to other off-diagonal entries
    memcpy(ptr2, ptr1, sizeof(Type)*dim2);

    //accumulate to diagonal entries
    ptrd1 = crs.A->GetPointer(n1, n1);
    ptrd2 = crs.A->GetPointer(n2, n2);
    for(i = 0; i < dim2; i++){
      ptrd1[i] -= ptr1[i];
      ptrd2[i] -= ptr1[i];
    }
  }
  
  //now do the boundary edges (only the parallel ones)
  for(eid = nbedge; eid < nbedge+ngedge; eid++){
    avec = m->bedges[eid].a;
    area = avec[3];
    n1 = m->bedges[eid].n[0];
    n2 = m->bedges[eid].n[1];
    if(!m->IsGhostNode(n2)){
      std::cerr << "WARNING: something wrong in lin. elastic" << std::endl;
      std::cerr << "Nnode: " << nnode << " Gnode: " << gnode 
		<< " Node1: " << n1 << " Node2: " << n2 
		<< " Edgeid: " << eid << std::endl;
      
      continue;
    }
    //TODO: generalize this so we have the volume assoc. with parallel nodes for movement
    //E = (aspectRatio[n1] + aspectRatio[n2])/2.0;
    E = aspectRatio[n1];

    //This section based on directional derivatives as given
    //Hyam's disseration p. 21
    ds2 = 0.0;
    dot = 0.0;
    for(j = 0; j < 3; j++){
      de[j] = m->xyz[n2*3 + j] - m->xyz[n1*3 + j];
      ds2 += de[j]*de[j];
    }
    for(j = 0; j < 3; j++){
      dot += de[j]*avec[j];
    }
    dot *= area;

    a11 = a22 = a33 = E*fact1;
    a12 = a13 = a21 = a23 = a31 = a32 = 0.5*E/nup1;
    t11 = t22 = t33 = E*fact2;
    t12 = t13 = t21 = t23 = t31 = t32 = 0.5*E/nup1;

    //off diagonal entry for row i1, column j1 -- dfQ(node1)/dQ(node2)
    //normal vector points from i1 to i2, differentiate edge flux w.r.t. to QR
    ptr1 = crs.A->GetPointer(n1, n2);
    //only real rows will have off-diagonal entries
    ptr1[0] = (a11+a12+a13)*dot/ds2;
    ptr1[3] = (t11*de[1]*avec[0] + t12*de[0]*avec[1])*area/ds2;
    ptr1[6] = (t11*de[2]*avec[0] + t13*de[0]*avec[2])*area/ds2;
    
    ptr1[1] = (t21*de[1]*avec[0] + t22*de[0]*avec[1])*area/ds2;
    ptr1[4] = (a21+a22+a23)*dot/ds2; 
    ptr1[7] = (t22*de[2]*avec[1] + t23*de[1]*avec[2])*area/ds2;
    
    ptr1[2] = (t31*de[2]*avec[0] + t33*de[0]*avec[2])*area/ds2;
    ptr1[5] = (t32*de[2]*avec[1] + t33*de[1]*avec[2])*area/ds2;
    ptr1[8] = (a31+a32+a33)*dot/ds2;	       

    //accumulate to diagonal entries
    ptrd1 = crs.A->GetPointer(n1, n1);
    for(i = 0; i < dim2; i++){
      ptrd1[i] -= ptr1[i];
    }
  }

  //set row to identity and zero rhs if node is on surface and locked
  SetBCLinearElastic(m, &crs, bc, dx);

  //setup rhs (equal to displacements)
  memcpy(crs.b, dx, (nnode*3)*sizeof(Type));

  //This seems to perform better with the relaxation SGS solver
  //GMRES does not seem to give accurate derivatives except with
  //huge # of search directions
#if 1
  //compute LU decomp for diagonals
  crs.A->PrepareSGS();
  //run SGS solver
  deltaDq = crs.SGS(iter, NULL, NULL, NULL, 1);
#else
  //run GMRES solver
  //use block diagonal preconditioning
  deltaDq = crs.GMRES(1, iter, 2, NULL, NULL, NULL);
#endif

  //copy back the displacements to the dx vector since we may use these to compute
  //the node velocities in a moving mesh situation
  memcpy(dx, crs.x, sizeof(Type)*(nnode*3));
  
  //update xyz coords locally
  for(i = 0; i < nnode*3; i++){
    m->xyz[i] += crs.x[i];
  }

  //update xyz coords across processes
  m->p->UpdateXYZ(m->xyz);

  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);
  std::cout << "Linear elastic solve - D||Ax-b||: "  << deltaDq << std::endl;

  delete [] aspectRatio;

  return;
}

template <class Type>
void SetBCLinearElastic(Mesh<Type>* m, CRS<Type>* crs, BoundaryConditions<Real>* bc, Type* dx)
{
  Int i;
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  Int* fixed = new Int[nnode+gnode];
  Type* ptr;

  //we are going to allow symmetry plane nodes to move
  Int* symNodes;
  Int nSymNodes = GetNodesOnSymmetryPlanes(m, bc, &symNodes);

  //change nodes on design plane but not directly specified to moving
  Int* designNodes;
  Int nDesignNodes = GetNodesOnMovingBC(m, bc, &designNodes);  

  //we are going to fix all nodes on the farfield BC
  Int* farNodes;
  Int nFarNodes = GetNodesOnBCType(m, bc, &farNodes, FarField);

  //init all nodes as moving
  for(i = 0; i < nnode+gnode; i++){
    fixed[i] = 0;
  }
  //change volume nodes to moving
  for(i = 0; i < m->nvnodes; i++){
    fixed[m->vnodes[i]] = 0;
  }
  //change symmetry nodes to moving
  for(i = 0; i < nSymNodes; i++){
    fixed[symNodes[i]] = 0;
  }
  //change nodes on the farfield BC to fixed - i.e. static
  for(i = 0; i < nFarNodes; i++){
    fixed[farNodes[i]] = 1;
  }
  //change nodes on a design flagged bc to fixed - i.e. prescribed
  for(i = 0; i < nDesignNodes; i++){
    fixed[designNodes[i]] = 1;
  }

  //change diag jacobian to the identity matrix for fixed nodes
  //and zero the off-diagonal row entries
  for(i = 0; i < nnode; i++){
    if(fixed[i]){
      crs->BlankMatrixRow(i);
      ptr = crs->A->GetPointer(i, i);
      ptr[0] = 1.0;
      ptr[4] = 1.0;
      ptr[8] = 1.0;
    }
  }

  delete [] fixed;
  delete [] symNodes;
  delete [] designNodes;

  return;
}

template <class Type>
void MoveMeshDesign(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc)
{
  //number of sweeps to move mesh by dx total
  Int nsweeps = 1;
  Int sweep, i;

#if 0
  //use scaled weighting
  Int weighting = 0;
  Real tol = 1.0e-11;
  Int max_iter = 20;
  Bool no_movement_dir[3];
  no_movement_dir[0] = no_movement_dir[1] = 0;
  //don't allow z-axis smoothing
  no_movement_dir[2] = 1;
  //smooth the mesh
  SmoothMeshLaplacian(&m, &bc, tol, max_iter, no_movement_dir, weighting);
#else
  Int iter = 1000;
  //we need this allocation to be big enough to store ghost node movement as well
  Int nnode = m->GetNumNodes();
  Type* dx = new Type[nnode*3];
  std::string dfname = casename;
  Int ndv = GetNdvDesignFile(dfname);
  Int beta;
  
  Type* dxyz = new Type[nnode*3];
  //initialize displacements to zero
  MemBlank(dx, nnode*3);
  
  for(beta = 0; beta < ndv; beta++){
    MemBlank(dxyz, nnode*3);
    Compute_dX_Surface(casename, m, bc, dxyz, beta); 
    for(i = 0; i < nnode*3; i++){
      dx[i] += dxyz[i];
    }
  }

  //divide total dx by nsweeps
  for(i = 0; i < nnode*3; i++){
    dx[i] /= (Type)nsweeps;
  }

  for(sweep = 0; sweep < nsweeps; sweep++){

    //move the mesh
    MoveMeshLinearElastic(m, bc, dx, iter);
    
    //recalculate metrics for next pass
    m->CalcMetrics();
  }
  
  delete [] dxyz;
  delete [] dx;
#endif
  
  return;
}
