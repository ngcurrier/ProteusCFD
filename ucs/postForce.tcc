#include "mesh.h"
#include "eqnset.h"
#include "bc.h"
#include "parallel.h"
#include "param.h"
#include "geometry.h"
#include "gradient.h"
#include "composite.h"

template <class Type>
PostForcePlugin<Type>::PostForcePlugin(SolutionSpace<Type>& space) :
  PostSurfaceIntegral<Type>(space, "Forces", Post_Force_Kernel),
  cf(NULL), cp(NULL)
{
  //allocate the correct number of composite bodies
  bodies = new CompositeBody<Type>[space.bc->num_bodies+1];

  Mesh<Type>* m = space.m;
  Int nbedge = m->GetNumBoundaryEdges();
  cf = new Type[nbedge];
  cp = new Type[nbedge];

  if(cp == NULL || cf == NULL || bodies == NULL){
    std::cerr << "Solution memory allocation error -- PostForcePlugin::PostForcePlugin()!!" 
	       << std::endl;
  }
}

template <class Type>
PostForcePlugin<Type>::~PostForcePlugin()
{
  //cleanup
  delete [] cp;
  delete [] cf;
  delete [] bodies; //bury the bodies...
}

template <class Type>
void PostForcePlugin<Type>::WriteTabularHeader()
{
  //todo...

}

template <class Type>
void PostForcePlugin<Type>::Report() const
{
  //todo...
}

template <class Type>
void Post_Force_Kernel(B_KERNEL_ARGS)
{
  Int i, j;

  EqnSet<Type>* eqnset = space->eqnset;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = eqnset->param;
  Mesh<Type>* m = space->m;
  PostForcePlugin<Type>* plugin = static_cast<PostForcePlugin<Type>*>(space->GetPlugin("Forces"));

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* Qsurf = &space->q[left_cv*nvars];
  Int bcType = bc->GetBCType(factag);
  Type gamma = param->gamma;
  Type cp = eqnset->GetCp(Qsurf, gamma);
  Type p = eqnset->GetPressure(Qsurf);
  if(m->IsGhostNode(right_cv)){
    //if we are on a parallel boundary don't compute anything here
    return;
  }
  plugin->cp[eid] = cp;
  
  //find number of terms passed in via the gradient
  std::vector<Int> gradLoc;
  Int nterms = eqnset->GetGradientsLocation(gradLoc);
  Type* grad = &space->qgrad[nterms*3*left_cv + 0];

  Type* tforces = (Type*)alloca(sizeof(Type)*3);
  Type* rmoment = (Type*)alloca(sizeof(Type)*3);
  Type* rpos = (Type*)alloca(sizeof(Type)*3);
  Type* momentCG;

  for(i = 1; i <= bc->GetNumBodies(); i++){
    if(plugin->bodies[i].SurfIsPart(factag)){
      //get the pt of the composite body to sum moments about
      momentCG = plugin->bodies[i].momentPt;
      //create position vector to the center of the face where the force acts
      Subtract(&m->cg[right_cv*3], momentCG, rpos);
    } 
    else{
      continue;
    }

    //add forces to temp vector of appropriate boundary object
    for(j = 0; j < 3; j++){
      tforces[j] = p*avec[j]*avec[3];
    }
    CrossProduct(rpos, tforces, rmoment);
    for(j = 0; j < 3; j++){
      plugin->bodies[i].forces[j] += tforces[j];
      plugin->bodies[i].moments[j] += rmoment[j];
    }

    //add the viscous forces
    if(bcType == Proteus_NoSlip && param->viscous){
      Type* stress = (Type*)alloca(sizeof(Type)*3);
      Type mu, rho, nu;
      mu = eqnset->ComputeViscosity(Qsurf);
      rho = eqnset->GetDensity(Qsurf);
      nu = mu/rho;
      eqnset->ComputeStressVector(grad, avec, mu, stress);
      for(j = 0; j < 3; j++){
	tforces[j] = stress[j]*avec[3];
      }
      CrossProduct(rpos, tforces, rmoment);
      for(j = 0; j < 3; j++){
	plugin->bodies[i].vforces[j] += tforces[j];
	plugin->bodies[i].vmoments[j] += rmoment[j];
      }
    }
  }
}
