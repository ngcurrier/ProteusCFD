#include "solutionSpace.h"
#include "alloca.h"
#include "exceptions.h"

template <class Type>
Int Driver(SolutionSpace<Type>* space, Kernel<Type> &kernel, Int scatterSize, 
	   void* custom)
{
  Int i, j, k;
  Int indx1, indx2;
  Int inode, eid;
  Int left_cv, right_cv;
  Mesh<Type>* m = space->m;
  
#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(space, m, custom, scatterSize) \
  private(i, j, k, inode, eid, indx1, indx2, left_cv, right_cv) 
#endif
  {
    Type* avec;
    Type* ptrL = NULL;
    Type* ptrR = NULL;
    Int n;
    Type* tempspaceL = new Type[scatterSize];
    Type* tempspaceR = new Type[scatterSize];
    Type* ve = (Type*)alloca(sizeof(Type)*3);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(i = 0; i < m->nedge; i++){
      eid = i;
      std::vector<pthread_mutex_t*> mutexes;
      left_cv = m->edges[eid].n[0];
      right_cv = m->edges[eid].n[1];
      avec = m->edges[eid].a;

      //get the velocity of the face
      for(k = 0; k < 3; k++){
	//take the average velocity at the edge face
	//TODO: generalize this for non-midpoint edges
	ve[k] = 0.5*(m->nv[left_cv*3 + k] + m->nv[right_cv*3 + k]);
      }
      //compute the velocity dotted with the edge normal, we compute the negative
      //b/c of the way it appears in the flux
      Type vdotn = -DotProduct(avec, ve);
      ptrL = ptrR = NULL;
      kernel(space, inode, left_cv, right_cv, avec, vdotn, &ptrL, &ptrR, 
	     tempspaceL, tempspaceR, &n, custom, eid);
      
      Type x = m->xyz[left_cv*3 + 0];
      Type y = m->xyz[left_cv*3 + 1];
      Type z = m->xyz[left_cv*3 + 2];
      std::stringstream ss;
      ss << "Coordinates for soft abort (kernel): x " << x << " y " << y << " z " << z;
      Abort.CheckForSoftAbort(ss.str());
      
      DriverScatter(left_cv, right_cv, ptrL, ptrR, tempspaceL, tempspaceR, n, mutexes);
    }
    delete [] tempspaceL;
    delete [] tempspaceR;
  }

  return (0);
}

template <class Type>
Int Bdriver(SolutionSpace<Type>* space, Kernel<Type> &bkernel, Int scatterSize, 
	    void* custom)
{
  Int eid;
  Int left_cv, right_cv;
  Type* avec;
  Type* ptrL = NULL;
  Type* ptrR = NULL;
  Int n, k;
  Int factag;
  Mesh<Type>* m = space->m;

#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(space, m, custom, scatterSize)\
  firstprivate(ptrL, ptrR)			   	       \
  private(eid, left_cv, right_cv, avec, n, factag, k) 
#endif
  {
    Type* tempspaceL = new Type[scatterSize];
    Type* tempspaceR = new Type[scatterSize];
    Type* ve = (Type*)alloca(sizeof(Type)*3);
    //loop over all the boundary edges in order
    //order not considered at this point
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(eid = 0; eid < m->nbedge+m->ngedge; eid++){
      std::vector<pthread_mutex_t*> mutexes;
      left_cv = m->bedges[eid].n[0];
      right_cv = m->bedges[eid].n[1];
      factag = m->bedges[eid].factag;
      avec = m->bedges[eid].a;
      //get the velocity of the face
      if(m->IsGhostNode(right_cv)){
	for(k = 0; k < 3; k++){
	  //take the average velocity at the edge face
	  //TODO: generalize this for non-midpoint edges
	  ve[k] = 0.5*(m->nv[left_cv*3 + k] + m->nv[right_cv*3 + k]);
	}
      }
      else{
	for(k = 0; k < 3; k++){
	  ve[k] = m->nv[left_cv*3 + k];
	}
      }
      //compute the velocity dotted with the edge normal, we compute the negative
      //b/c of the way it appears in the flux
      Type vdotn = -DotProduct(avec, ve);
      ptrL = ptrR = NULL;
      bkernel(space, left_cv, left_cv, right_cv, avec, vdotn, ve, &ptrL, &ptrR, 
	      tempspaceL, tempspaceR, &n, custom, eid, factag);
      
      Type x = m->xyz[left_cv*3 + 0];
      Type y = m->xyz[left_cv*3 + 1];
      Type z = m->xyz[left_cv*3 + 2];
      std::stringstream ss;
      ss << "Coordinates for soft abort(bkernel): x " << x << " y " << y << " z " << z;
      Abort.CheckForSoftAbort(ss.str());

      DriverScatter(left_cv, right_cv, ptrL, ptrR, tempspaceL, tempspaceR, n, mutexes);
    }
    delete [] tempspaceL;
    delete [] tempspaceR;
  }
  return (0);
}

//This is exclusively used to set boundary conditions since 
//driver scattering makes no sense in that context
template <class Type>
Int BdriverNoScatter(SolutionSpace<Type>* space, Kernel<Type> &bkernel, Int scatterSize, 
		     void* custom)
{
  Int eid;
  Int left_cv, right_cv;
  Type* avec;
  Int n, k;
  Int factag;
  Mesh<Type>* m = space->m;

#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(space, m, custom, scatterSize) \
  private(eid, left_cv, right_cv, avec, n, factag, k) 
#endif
  {
    Type* ptrL = NULL;
    Type* ptrR = NULL;
    Type* tempspaceL = new Type[scatterSize];
    Type* tempspaceR = new Type[scatterSize];
    Type* ve = (Type*)alloca(sizeof(Type)*3);
    //loop over all the boundary edges in order
    //order not considered at this point#ifdef _OPENMP
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(eid = 0; eid < m->nbedge+m->ngedge; eid++){
      left_cv = m->bedges[eid].n[0];
      right_cv = m->bedges[eid].n[1];
      factag = m->bedges[eid].factag;
      avec = m->bedges[eid].a;
      //get the velocity of the face
      if(m->IsGhostNode(right_cv)){
	for(k = 0; k < 3; k++){
	  //take the average velocity at the edge face
	  //TODO: generalize this for non-midpoint edges
	  ve[k] = 0.5*(m->nv[left_cv*3 + k] + m->nv[right_cv*3 + k]);
	}
      }
      else{
	for(k = 0; k < 3; k++){
	  ve[k] = m->nv[left_cv*3 + k];
	}
      }
      //compute the velocity dotted with the edge normal, we compute the negative
      //b/c of the way it appears in the flux
      Type vdotn = -DotProduct(avec, ve);
      ptrL = ptrR = NULL;
      bkernel(space, left_cv, left_cv, right_cv, avec, vdotn, ve, &ptrL, &ptrR, 
	      tempspaceL, tempspaceR, &n, custom, eid, factag);

      Type x = m->xyz[left_cv*3 + 0];
      Type y = m->xyz[left_cv*3 + 1];
      Type z = m->xyz[left_cv*3 + 2];
      std::stringstream ss;
      ss << "Coordinates for soft abort(bkernel - no scatter): x " << x << " y " << y << " z " << z;
      Abort.CheckForSoftAbort(ss.str());
    }
    delete [] tempspaceL;
    delete [] tempspaceR;
  }
  return (0);
}

template <class Type>
Int DriverNoScatter(SolutionSpace<Type>* space, Kernel<Type> &kernel, 
		    Int scatterSize, void* custom)
{
  Int i, j, k;
  Int indx1, indx2;
  Int inode, eid;
  Int left_cv, right_cv;
  Mesh<Type>* m = space->m;

#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(none) shared(space, m, custom, scatterSize) \
  private(i, j, k, inode, eid, indx1, indx2, left_cv, right_cv) 
#endif
  {
    Type* avec;
    Type* ptrL = NULL;
    Type* ptrR = NULL;
    Int n;
    Type* tempspaceL = new Type[scatterSize];
    Type* tempspaceR = new Type[scatterSize];
    Type* ve = (Type*)alloca(sizeof(Type)*3);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(i = 0; i < m->nedge; i++){
      eid = i;
      left_cv = m->edges[eid].n[0];
      right_cv = m->edges[eid].n[1];
      avec = m->edges[eid].a;
      //get the velocity of the face
      for(k = 0; k < 3; k++){
	//take the average velocity at the edge face
	//TODO: generalize this for non-midpoint edges
	ve[k] = 0.5*(m->nv[left_cv*3 + k] + m->nv[right_cv*3 + k]);
      }
      //compute the velocity dotted with the edge normal, we compute the negative
      //b/c of the way it appears in the flux
      Type vdotn = -DotProduct(avec, ve);
      ptrL = ptrR = NULL;
      kernel(space, inode, left_cv, right_cv, avec, vdotn, &ptrL, &ptrR, 
	     tempspaceL, tempspaceR, &n, custom, eid);
      Type x = m->xyz[left_cv*3 + 0];
      Type y = m->xyz[left_cv*3 + 1];
      Type z = m->xyz[left_cv*3 + 2];
      std::stringstream ss;
      ss << "Coordinates for soft abort (kernel - no scatter): x " << x << " y " << y << " z " << z;
      Abort.CheckForSoftAbort(ss.str());
    }
    delete [] tempspaceL;
    delete [] tempspaceR;
  }

  return (0);
}


//used for scattering on driver side with openmp
template <class Type>
void DriverScatter(Int left_cv, Int right_cv, Type* ptrL, Type* ptrR, Type* valL, Type* valR, 
		   Int n, std::vector<pthread_mutex_t*>& mutexes)
{
  Int i;
#ifdef _OPENMP
  LockVars(mutexes);
  if(ptrR != NULL){
    for(i = 0; i < n; i++){
      ptrR[i] += valR[i];
    }
  }
  if(ptrL != NULL){
    for(i = 0; i < n; i++){
      ptrL[i] += valL[i];
    }
  }
  UnlockVars(mutexes);
#else

  if(ptrR != NULL){
    for(i = 0; i < n; i++){
      ptrR[i] += valR[i];
    }
  }
  if(ptrL != NULL){
    for(i = 0; i < n; i++){
      ptrL[i] += valL[i];
    }
  }
#endif

}
