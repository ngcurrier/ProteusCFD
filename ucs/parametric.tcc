#include "bc.h"
#include "macros.h"
#include "mesh.h"
#include "portFileio.h"
#include "exceptions.h"
#include "ffd.h"
#include <sstream>
#include <cmath>

template <class Type>
void Compute_dX_Surface(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			Type* dxyz, Int beta)
{
  Int i;
  Int parType;
  std::string designName = casename;
  Int nnode = m->GetNumNodes();

  Int ndv = GetNdvDesignFile(designName);
  Real* x = new Real[ndv];
  Real* bounds = new Real[ndv*2];
  Real f;
  Real* grad = new Real[ndv];
  Int* dvType = new Int[ndv];
  Int pt;

  MemBlank(dxyz, nnode*3);

  ReadDesignFile(designName, &ndv, x, bounds, &f, grad, dvType); 

  parType = dvType[beta];
  
  //direct point movement
  if(parType == 0 || parType == 1 || parType == 2){
    
    Abort << "Direct node movement not function b/c nodeorigtolocal calls no longer possible";

    Int localNode;
    //std::string pointFile = PT_FILE;
    //Int npts = GetNptsPointFile(pointFile);
    Int npts = 0;
    Int* pts = new Int[npts];
    Real* xyz = new Real[npts*3];
    Real* ptdxyz = new Real[npts*3];
    Int* betaid = new Int[npts];

    //ReadPointMoveFile(pointFile, npts, ndv, pts, xyz, ptdxyz, betaid);
    for(i = 0; i < npts; i++){
      if(betaid[i] == beta){
	pt = pts[i];
		
	//localNode = m->NodeOrigToLocal(pt);

	//make sure the node is on this process
	if(localNode >= 0){
	  if(parType == 0){
	    dxyz[localNode*3 + 0] = x[beta];
	  }
	  else if(parType == 1){
	    dxyz[localNode*3 + 1] = x[beta];
	  }
	  else if(parType == 2){
	    dxyz[localNode*3 + 2] = x[beta];
	  }
	  else{
	    //should never get here...
	  }
	}
      }
    }
    
    delete [] pts;
    delete [] xyz;
    delete [] ptdxyz;
    delete [] betaid;
  }
  //hicks-henne function
  //assume these are placed in design file in order
  else if(parType == 3){
    //if(parType == 3 || parType == 4 || parType == 5){
    Int* pts;
    Int pt;
    Int spts = GetNodesOnMovingBC(m, bc, beta, &pts);
    Type* sxyz = new Type[spts*3];
    Type* sx = new Type[spts];
    Type* sy = new Type[spts];

    //get coordinates of points on surface
    m->GetCoordsForPoints(sxyz, pts, spts);
    //strip off x coordinates for the points on the surface
    for(i = 0; i < spts; i++){
      sx[i] = sxyz[i*3 + 0];
    }
    //find largest x value among all processes for normalization
    PObj<Type>* p = m->p;
    Int np = p->GetNp();
    Int rank = p->GetRank();
    Type* max = (Type*)alloca(sizeof(Type)*np);
    for(i = 0; i < np; i++){
      max[i] = 0.0;
    }
    for(i = 0; i < spts; i++){
      max[rank] = MAX(max[rank], sx[i]);
    }
    //reduce across all processes, MPI_MAX is not valid for complex types
    //so we have to work around it... very not cool
    MPI_Datatype mpit;
    mpit = MPI_GetType(max[0]);
    MPI_Allgather(MPI_IN_PLACE, 1, mpit, max, 1, mpit, MPI_COMM_WORLD); 

    for(i = 0; i < np; i++){
      max[0] = MAX(max[i], max[0]);
    }

    //normalize x coords - HH function requires it
    NormalizeX(sx, spts, max[0]);

    //a - magnitude (3)
    if(parType == 3){
      Hicks_Henne((Type)x[beta], (Type)x[beta+1], (Type)x[beta+2], sx, sy, spts);
    }
    //t1 - max point of bump (4) .. if plotted on 0-1.0 should be greater than ~1.1 - ~1000.0, with 
    //                               smaller values place peak near trailing edge, near x = 1.0
    else if(parType == 4){
      Hicks_Henne((Type)x[beta-1], (Type)x[beta], (Type)x[beta+1], sx, sy, spts);
    }
    //t2 - width of bump (5) .. if plotted on 0.0-1.0 should be 0 - ~2000 , larger is narrower
    else if(parType == 5){
      Hicks_Henne((Type)x[beta-2], (Type)x[beta-1], (Type)x[beta], sx, sy, spts);
    }
    else{
      //should never get here...
    }

    //load up displacement array
    for(i = 0; i < spts; i++){
      pt = pts[i];
      dxyz[pt*3 + 1] = sy[i];
    }
    delete [] pts;
    delete [] sxyz;
    delete [] sx;
    delete [] sy;
  }
  //this is the fall through from above... do nothing
  else if(parType == 4 || parType == 5){

  }
  //this is the x, y, and z motion of a whole boundary condition
  else if(parType == 6 || parType == 7 || parType == 8){
    Int* pts;
    Int spts = GetNodesOnMovingBC(m, bc, beta, &pts);

    Type* sxyz = new Type[spts*3];

    Type dx, dy, dz;
    dx = dy = dz = 0.0;

    if(parType == 6){
      dx = x[beta];
    }
    else if(parType == 7){
      dy = x[beta];
    }
    else if(parType == 8){
      dz = x[beta];
    }
    else{
      //should never get here...
    }

    //get coordinates of points on surface
    m->GetCoordsForPoints(sxyz, pts, spts);

    //load up displacement array
    for(i = 0; i < spts; i++){
      pt = pts[i];
      dxyz[pt*3 + 0] = dx;
      dxyz[pt*3 + 1] = dy;
      dxyz[pt*3 + 2] = dz;
    }

    delete [] pts;
    delete [] sxyz;
  }
  //this is the gaussian plume location parameters
  else if(parType == 9 || parType == 10){
    //do nothing
  }
  else{
    //should never be here
    std::stringstream ss;
    ss << "In Compute_dX_Surface() and design parameter type ";
    ss << parType;
    ss << " not found!";
    
    Abort << ss.str();
  }

  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;

  return;
}

template <class Type>
void Compute_dXdB_Surface(std::string casename, Mesh<Type>* m, BoundaryConditions<Real>* bc, 
			  Type* dxdbSurf, Int beta)
{
  Int i;
  
  std::string designName = casename;
  
  Int ndv = GetNdvDesignFile(designName);
  Type* x = new Type[ndv];
  Type* bounds = new Type[ndv*2];
  Type f;
  Type* grad = new Type[ndv];
  Int* dvType = new Int[ndv];
  Int parType;
  Int nnode = m->GetNumNodes();

  MemBlank(dxdbSurf, nnode*3);

  ReadDesignFile(designName, &ndv, x, bounds, &f, grad, dvType); 
  
  parType = dvType[beta];

  if(parType == 0 || parType == 1 || parType == 2){
    Int localNode;
    //std::string pointFile = PT_FILE;
    //Int npts = GetNptsPointFile(pointFile);
    Int npts = 0;
    Int* pts = new Int[npts];
    Type* xyz = new Type[npts*3];
    Type* ptdxyz = new Type[npts*3];
    Int* betaid = new Int[npts];

    //ReadPointMoveFile(pointFile, npts, ndv, pts, xyz, ptdxyz, betaid);
    for(i = 0; i < npts; i++){
      if(betaid[i] == beta){
	localNode = pts[i];
	//make sure the node is on this process
	if(localNode >= 0){
	  if(parType == 0){
	    dxdbSurf[localNode*3 + 0] = 1.0;
	  }
	  else if(parType == 1){
	    dxdbSurf[localNode*3 + 1] = 1.0;
	  }
	  else if(parType == 2){
	    dxdbSurf[localNode*3 + 2] = 1.0;
	  }
	  else{
	    //should never get here...
	  }
	}
      }
    }
    delete [] pts;
    delete [] xyz;
    delete [] ptdxyz;
    delete [] betaid;
  }
  //hicks_henne functions
  else if(parType == 3 || parType == 4 || parType == 5){
    Int* pts;
    Int pt;
    Int spts = GetNodesOnMovingBC(m, bc, beta, &pts);
    Type* sxyz = new Type[spts*3];
    Type* sx = new Type[spts];
    Type* sy = new Type[spts];
    
    //strip off x coordinates of points on surface
    m->GetCoordsForPoints(sxyz, pts, spts);
    for(i = 0; i < spts; i++){
      sx[i] = sxyz[i*3 + 0];
    }
    //find largest x value among all processes for normalization
    Type max = 0.0;
    for(i = 0; i < spts; i++){
      max = MAX(max, sx[i]);
    }
    //reduce across all processes
    MPI_Datatype mpit;
    mpit = MPI_GetType(max);
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, mpit, MPI_MAX, MPI_COMM_WORLD);

    //normalize x coords - HH function requires it
    NormalizeX(sx, spts, max);

    //a - magnitude (3)
    if(parType == 3){
      Da_Hicks_Henne(x[beta], x[beta+1], x[beta+2], sx, sy, spts);
    }
    //t1 - max point of bump (4)
    else if(parType == 4){
      Dt1_Hicks_Henne(x[beta-1], x[beta], x[beta+1], sx, sy, spts);
    }
    //t2 - width of bump (5)
    else if(parType == 5){
      Dt2_Hicks_Henne(x[beta-2], x[beta-1], x[beta], sx, sy, spts);
    }
    else{
      //should never get here...
    }
    
    //load up displacement array
    for(i = 0; i < spts; i++){
      pt = pts[i];
      dxdbSurf[pt*3 + 1] = sy[i];
    }
    delete [] pts;
    delete [] sxyz;
    delete [] sx;
    delete [] sy;
  }
  //this is the x, y, and z motion of a whole boundary condition
  else if(parType == 6 || parType == 7 || parType == 8){
    Int pt;
    Int* pts;
    Int spts = GetNodesOnMovingBC(m, bc, beta, &pts);

    Type* sxyz = new Type[spts*3];

    Type dx, dy, dz;
    dx = dy = dz = 0.0;

    if(parType == 6){
      dx = 1.0;
    }
    else if(parType == 7){
      dy = 1.0;
    }
    else if(parType == 8){
      dz = 1.0;
    }
    else{
      //should never get here...
    }

    //get coordinates of points on surface
    m->GetCoordsForPoints(sxyz, pts, spts);

    //load up displacement array
    for(i = 0; i < spts; i++){
      pt = pts[i];
      dxdbSurf[pt*3 + 0] = dx;
      dxdbSurf[pt*3 + 1] = dy;
      dxdbSurf[pt*3 + 2] = dz;
    }

    delete [] pts;
    delete [] sxyz;
  }
  //this is the x and y location for the gaussian plume
  else if(parType == 9 || parType == 10){
    //no movement
  }
  else{
    //should never be here
    std::stringstream ss;
    ss << "In Compute_dXdB_Surface() and design parameter type ";
    ss << parType;
    ss << " not found!";
    
    Abort << ss.str();
  }

  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;
   
  return;
}

template <class Type>
void Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n)
{
  Int i;
  Type c1 = log(5.0)/log(t1);
  Type pi = acos(-1.0);

  for(i = 0; i < n; i++){
    y[i] = a*pow(sin(pi*pow(x[i],c1)), t2);
  }

  return;
}

template <class Type>
void Da_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n)
{
  Int i;
  Type c1 = log(5)/log(t1);
  Type pi = acos(-1.0);

  for(i = 0; i < n; i++){
    //derivative from maple
    y[i] = pow(sin(pi*pow(x[i],c1)), t2);
  }

  return;
}

template <class Type>
void Dt1_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n)
{
  Int i;
  Type c1 = log(5)/log(t1);
  Type P;
  Type pi = acos(-1.0);

  for(i = 0; i < n; i++){
    P = pow(x[i],c1);
    //derivative from maple
    y[i] = (-a*pow(sin(pi*P),t2-1.0)*t2*cos(pi*P)*pi*P*log(5.0)*log(x[i]))/
    (log(t1*t1)*t1);
  }

  return;
}

template <class Type>
void Dt2_Hicks_Henne(Type a, Type t1, Type t2, Type* x, Type* y, Int n)
{
  Int i;
  Type c1 = log(5)/log(t1);
  Type P;
  Type pi = acos(-1.0);

  for(i = 0; i < n; i++){
    P = pow(x[i],c1);
    //derivative from maple
    y[i] = a*pow(sin(pi*P), t2)*log(sin(pi*P));
  }

  return;
}


template <class Type>
Type NormalizeX(Type* x, Int n, Type max)
{
  Int i;

  if(max == 0.0){
    //find max point
    max = x[0];
    for(i = 0; i < n; i++){
      max = MAX(x[i], max);
    }
  }

  //normalize points
  for(i = 0; i < n; i++){
    x[i] /= max;
  }
  return max;
}
