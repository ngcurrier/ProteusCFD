#include "driver.h"
#include "mesh.h"
#include "eqnset.h"
#include "bc.h"
#include "parallel.h"
#include "param.h"
#include "geometry.h"
#include "gradient.h"
#include "solutionSpace.h"
#include "composite.h"

template <class Type>
Forces<Type>::Forces(SolutionSpace<Type>* space)
{
  //set all appropriate pointers
  this->space = space;

  Param<Type>* param = space->param;

  //open the forces outfile for writing
  std::string forcesFileName = param->path+param->spacename + ".cl";
  if(space->p->GetRank() == 0){
    //if not using restart or using restart but not preserving counters, open a new file
    if(!param->useRestart || !param->preserveRestartCounters){
      fout.open(forcesFileName.c_str());
    }
    //if restarting, append to the current file
    else{
      fout.open(forcesFileName.c_str(), std::ios::app);
    }
  }
  //set output print formats
  fout.setf(std::ios::scientific);
  fout.precision(6);

  forcesAlloc = false;

  //we need to push these back on the scalar fields for typical section analysis
  space->AddScalarField("CL");
  space->AddScalarField("CM");
}

template <class Type>
Forces<Type>::~Forces()
{
  fout.close();
  KillForcesMemory();

}


template <class Type>
void Forces<Type>::Report()
{
  ReportCl();
  ReportYpCf();
}

template <class Type>
void Forces<Type>::ReportCl()
{
  Int i;
  PObj<Type>* p = space->p;
  Param<Type>* param = space->param;

  Int iter = space->iter;
  Type liftdir[3];
  Type dragdir[3];
  Type lift;
  Type drag;
  Type moment;

  for(i = 0; i < 3; i++){
    liftdir[i] = param->liftdir[i];
    dragdir[i] = param->dragdir[i];
  }

  Type* momdir;
  if(p->GetRank() == 0){
    fout << iter << ": " << std::endl;
    for(i = 1; i <= num_bodies; i++){
      momdir = bodies[i].momentAxis;
      std::string name = bodies[i].name;
      lift = DotProduct(liftdir, bodies[i].forces);
      drag = DotProduct(dragdir, bodies[i].forces);
      moment = DotProduct(momdir, bodies[i].moments);
      fout << "Force vector for body [" << i << " - " << name << "]: " 
	   << bodies[i].forces[0] << " " << bodies[i].forces[1] << " " 
	   << bodies[i].forces[2] << std::endl;
      fout << "Viscous force vector for body [" << i << " - " << name << "]: " 
	   << bodies[i].vforces[0] << " " << bodies[i].vforces[1] << " " 
	   << bodies[i].vforces[2] << std::endl;
      fout << "Moment vector for body[" << i << " - " << name << "]: " 
	   << bodies[i].moments[0] << " " << bodies[i].moments[1] << " " 
	   << bodies[i].moments[2] << std::endl;
      fout << "Viscous moment vector for body [" << i << " - " << name << "]: " 
	   << bodies[i].vmoments[0] << " " << bodies[i].vmoments[1] << " " 
	   << bodies[i].vmoments[2] << std::endl;
      fout << "Lift force for body[" << i << " - " << name << "]: " 
	   << lift << std::endl;   
      fout << "Moment for body[" << i << " - " << name << "]: " 
	   << moment << std::endl;
      fout << "Lift coefficient for body[" << i << " - " << name << "]: " 
	   << bodies[i].cl << std::endl;
      fout << "Moment coefficient for body[" << i << " - " << name << "]: " 
	   << bodies[i].cm << std::endl;
      fout << "Drag force for body[" << i << " - " << name << "]: " 
	   << drag << std::endl;
      fout << "Drag coefficient for body[" << i << " - " << name << "]: " 
	   << bodies[i].cd << std::endl;
    } 
  }
}


//compute total force acting on surface 
template <class Type>
void FORCE_Kernel(B_KERNEL_ARGS)
{
  Int i, j;

  EqnSet<Type>* eqnset = space->eqnset;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = eqnset->param;
  Mesh<Type>* m = space->m;
  Forces<Type>* forces = space->forces;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* Qsurf = &space->q[left_cv*nvars];
  Int bcId = bc->GetBCId(factag);
  Type gamma = param->gamma;
  Type cp = eqnset->GetCp(Qsurf, gamma);
  Type p = eqnset->GetPressure(Qsurf);
  if(m->IsGhostNode(right_cv)){
    //if we are on a parallel boundary don't compute anything here
    return;
  }
  forces->cp[eid] = cp;
  
  //find number of terms passed in via the gradient
  std::vector<Int> gradLoc;
  Int nterms = eqnset->GetGradientsLocation(gradLoc);
  Type* grad = &space->qgrad[nterms*3*left_cv + 0];

  Type* tforces = (Type*)alloca(sizeof(Type)*3);
  Type* rmoment = (Type*)alloca(sizeof(Type)*3);
  Type* rpos = (Type*)alloca(sizeof(Type)*3);
  Type* momentCG;

  for(i = 1; i <= forces->num_bodies; i++){
    if(forces->bodies[i].SurfIsPart(factag)){
      //get the pt of the composite body to sum moments about
      momentCG = forces->bodies[i].momentPt;
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
      forces->bodies[i].forces[j] += tforces[j];
      forces->bodies[i].moments[j] += rmoment[j];
    }

    //add the viscous forces
    if(bcId == NoSlip && param->viscous){
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
	forces->bodies[i].vforces[j] += tforces[j];
	forces->bodies[i].vmoments[j] += rmoment[j];
      }
    }
  }
}

template <class Type>
void ComputeSurfaceAreas(SolutionSpace<Type>* space, Int verbosity)
{
  Int i;
  BoundaryConditions<Real>* bc = space->bc;
  Forces<Type>* forces = space->forces;
  BCObj<Real>* bcs = bc->bcs;
  MPI_Datatype mpit;

  //get mpi datatype to send
  mpit = MPI_GetType(forces->surfArea[0]);


  //zero out mesh surface area for accumulation
  //recall here that bc_obj 0 is saved as a parallel (non-bc) condition
  for(i = 0; i <= forces->num_bcs; i++){
    forces->surfArea[i*3 + 0] = 0.0;
    forces->surfArea[i*3 + 1] = 0.0;
    forces->surfArea[i*3 + 2] = 0.0;
  } 
  //zero areas for the composite bodies
  for(i = 1; i <= forces->num_bodies; i++){
    forces->bodies[i].surfArea[0] = forces->bodies[i].surfArea[1] = forces->bodies[i].surfArea[2] = 0.0;
  }

  Kernel<Type> Surface_Area(SURFACE_AREA_Kernel);

  //call driver to loop over boundaries for accumulation
  BdriverNoScatter(space, Surface_Area, 3, NULL);

  //sum up areas over all processes
  MPI_Allreduce(MPI_IN_PLACE, forces->surfArea, 3*(forces->num_bcs+1), mpit, MPI_SUM, MPI_COMM_WORLD);

  if(verbosity > 0){
    //report surface areas
    std::cout << "Surface areas for identified boundaries:" << std::endl;
    std::cout << "========================================" << std::endl;
    for(i = 1; i <= forces->num_bcs; i++){
      //calculate area magnitude
      Type amag = Magnitude(&forces->surfArea[i*3]);
      std::cout << "Area for factag[" << bcs[i].factag << "] : " << amag << std::endl;   
    }
    std::cout << std::endl;
  }

  for(i = 1; i <= forces->num_bodies; i++){
    //sum up areas over all processes
    MPI_Allreduce(MPI_IN_PLACE, forces->bodies[i].surfArea, 3, mpit, MPI_SUM, MPI_COMM_WORLD);
  }

  if(verbosity > 0){
    //report surface areas
    std::cout << "Surface areas for composite bodies:" << std::endl;
    std::cout << "====================================" << std::endl;
    for(i = 1; i <= forces->num_bodies; i++){
      //calculate area magnitude
      Type amag = Magnitude(forces->bodies[i].surfArea);
      std::cout << "Area for body[" << i << "] : " << amag << std::endl;   
    }
    std::cout << std::endl;
  }
}

template <class Type>
void SURFACE_AREA_Kernel(B_KERNEL_ARGS)
{
  Int i, j;
  Forces<Type>* forces = space->forces;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = space->param;
  Int bcId = bc->bc_map[factag];
  Type dot;
  Type* liftdir = (Type*)alloca(sizeof(Type)*3);
  //if parallel boundary don't do anything.. 
  if(factag == 0){
    return;
  }

  forces->surfArea[bcId*3 + 0] += (CAbs(avec[0]*avec[3]));
  forces->surfArea[bcId*3 + 1] += (CAbs(avec[1]*avec[3]));
  forces->surfArea[bcId*3 + 2] += (CAbs(avec[2]*avec[3]));

  for(i = 0; i < 3; i++){
    liftdir[i] = param->liftdir[i];
  }

  for(i = 1; i <= forces->num_bodies; i++){
    if(forces->bodies[i].SurfIsPart(factag)){
#if 1
      dot = DotProduct(liftdir, avec);
      //we are only interested in the projected planform area
      //and only in the direction of the lift/drag forces
      if(real(dot) >= 0.0){
	for(j = 0; j < 3; j++){
	  forces->bodies[i].surfArea[j] += (CAbs(dot*avec[j]*avec[3]));
	}  
      }
#else
      for(j = 0; j < 3; j++){
	forces->bodies[i].surfArea[j] += (CAbs(avec[j]*avec[3]));
      }
#endif
    } 
  }
  
  return;
}

template <class Type>
void Forces<Type>::Compute()
{
  ComputeYpCf();
  ComputeCl();

  //set fields for typical section analysis
  //assume that the first body defined is what we are interested in
  space->GetScalarField("CL").SetField(bodies[1].cl);
  space->GetScalarField("CM").SetField(bodies[1].cm);
}

template <class Type>
void Forces<Type>::ComputeCl()
{
  Int i;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;
  Type* Qref = eqnset->Qinf;
  Type liftdir[3];
  Type dragdir[3];
  Type lift;
  Type drag;
  Type moment;

  for(i = 0; i < 3; i++){
    liftdir[i] = param->liftdir[i];
    dragdir[i] = param->dragdir[i];
  }

  Type* momdir;
  Type V = param->GetVelocity(space->iter);
  Type Mach = V;
  Type rho = eqnset->GetDensity(Qref);
  Type* area;
  Type amag;

  //Cl = L /(0.5 * rho * V^2 * A)
  //Cd = D /(0.5 * rho * V^2 * A)
  //Cm = M /(0.5 * rho * V^2 * A * Lref)
  for(i = 1; i <= num_bodies; i++){
    //get the pt of the composite body to sum moments about
    momdir = bodies[i].momentAxis;
    //compute the pressure lift/drag
    lift = DotProduct(liftdir, bodies[i].forces);
    drag = DotProduct(dragdir, bodies[i].forces);
    moment = DotProduct(momdir, bodies[i].moments);
    //compute the viscous lift/drag
    lift += DotProduct(liftdir, bodies[i].vforces);
    drag += DotProduct(dragdir, bodies[i].vforces);
    moment += DotProduct(momdir, bodies[i].vmoments);
    area = bodies[i].surfArea;
    amag = Magnitude(area);
    bodies[i].cl = lift / (0.5 * rho * Mach*Mach * amag);
    bodies[i].cd = drag / (0.5 * rho * Mach*Mach * amag);
    //assume that the ref_length is unity in the grid, otherwise, WRONG!
    //positive cm pitches the airfoil nose up
    bodies[i].cm = -moment / (0.5 * rho * Mach*Mach * amag * 1.0);
  } 
}

template <class Type>
void Forces<Type>::ComputeYpCf()
{
  Int i;
  MPI_Datatype mpit;
  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  Param<Type>* param = space->param;
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;

  //get mpi datatype to send
  Type tmp = 0.0;
  mpit = MPI_GetType(tmp);

  //blank the composite bodies forces, moments, etc.
  //remember bodies[0] is not used, generally
  for(i = 0; i <= num_bodies; i++){
    MemBlank(bodies[i].forces, 3);
    MemBlank(bodies[i].vforces, 3);
    MemBlank(bodies[i].moments, 3);
    MemBlank(bodies[i].vmoments, 3);
  }

  //zero cf, yp, qdot
  for(Int eid = 0; eid < m->GetNumBoundaryEdges(); eid++){
    yp[eid] = cf[eid] = qdot[eid] = 0.0;
  }

  //call driver to loop over boundaries for accumulation
  Kernel<Type> FORCE(FORCE_Kernel);
  BdriverNoScatter(space, FORCE, 3, NULL);

  //sum up forces, moments over all processes
  for(i = 1; i <= num_bodies; i++){
    MPI_Allreduce(MPI_IN_PLACE, bodies[i].forces, 3, mpit, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, bodies[i].vforces, 3, mpit, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, bodies[i].moments, 3, mpit, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, bodies[i].vmoments, 3, mpit, MPI_SUM, MPI_COMM_WORLD);
  }

  //call driver to loop over boundaries for accumulation of yplus, cf, qdot, etc.
  Kernel<Type> YpCf(YpCf_Kernel);
  BdriverNoScatter(space, YpCf, 0, NULL);
}

template <class Type>
void YpCf_Kernel(B_KERNEL_ARGS)
{
  Int i;

  Mesh<Type>* m = space->m;
  EqnSet<Type>* eqnset = space->eqnset;
  BoundaryConditions<Real>* bc = space->bc;
  Param<Type>* param = space->param;
  Forces<Type>* forces = space->forces;

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Type* Qsurf = &space->q[left_cv*nvars];
  Int bcId = bc->GetBCId(factag);

  Type* stress = (Type*)alloca(sizeof(Type)*3);
  Type tauw, mu, nu, Tsurf, rho;

  //find number of terms passed in via the gradient
  std::vector<Int> gradLoc;
  Int nterms = eqnset->GetGradientsLocation(gradLoc);

  //extract the correct gradient which has been precomputed
  Type* grad = &space->qgrad[nterms*3*left_cv + 0];
  Type* wallx = &m->xyz[left_cv*3];
  Type* dx = (Type*)alloca(sizeof(Type)*3);

  if(m->IsGhostNode(right_cv)){
    //we don't compute cf, etc. for parallel boundary edges
    //move along...
    return;
  }

  //only compute y+ and cf for viscous surfaces
  if(bcId == NoSlip && param->viscous){
    //we need to find the distance to the first grid point off the wall in the 
    //normal direction, can include parallel nodes
    Int indx, indx1, indx2;
    Type d = 0.0;
    indx1 = m->ipsp[left_cv];
    indx2 = m->ipsp[left_cv+1];
    Type dotmax = 0.0;
    for(indx = indx1; indx < indx2; indx++){
      Int pt = m->psp[indx];
      Type* ptx = &m->xyz[pt*3];
      Subtract(ptx, wallx, dx);
      Normalize(dx, dx);
      Type dot = -DotProduct(dx, avec);
      //find the most normal point off the wall, use that distance for yplus
      if(real(dot) >= real(dotmax)){
	d = Distance(ptx, wallx);
	dotmax = dot;
      }
    }
    mu = eqnset->ComputeViscosity(Qsurf);
    rho = eqnset->GetDensity(Qsurf);
    nu = mu/rho;
    eqnset->ComputeStressVector(grad, avec, mu, stress);
    tauw = Magnitude(stress);
    //multiply by ReTilde since tauw has that term squared in it
    //and yplus is not non-dimensionalized the same as the stress vector is
    forces->yp[eid] = d*sqrt(tauw/rho)/nu*eqnset->GetRe();
    forces->cf[eid] = eqnset->GetCf(tauw, rho);
  }
}

template <class Type>
void Forces<Type>::ReportYpCf()
{
  Int eid;
  Mesh<Type>* m = space->m;
  PObj<Type>* p = space->p;
  BoundaryConditions<Real>* bc = space->bc;
  Forces<Type>* forces = space->forces;

  Real maxYp = 0.0;
  Real minYp = 99999.99;
  Real maxYpLoc[3];
  Real minYpLoc[3];
  for(eid = 0; eid < m->GetNumBoundaryEdges(); eid++){
    Int factag = m->bedges[eid].factag;
    Int bcId = bc->GetBCId(factag);
    //we only care about the min and max yplus on viscous boundaries
    if(bcId == NoSlip){
      maxYp = MAX(maxYp, real(forces->yp[eid]));
      if(real(forces->yp[eid]) >= maxYp){
	maxYp = real(forces->yp[eid]);
	Int left_cv = m->bedges[eid].n[0];
	maxYpLoc[0] = real(m->xyz[left_cv*3 + 0]);
	maxYpLoc[1] = real(m->xyz[left_cv*3 + 1]);
	maxYpLoc[2] = real(m->xyz[left_cv*3 + 2]);
      }
      if(real(forces->yp[eid]) <= minYp){
	minYp = real(forces->yp[eid]);
	Int left_cv = m->bedges[eid].n[0];
	minYpLoc[0] = real(m->xyz[left_cv*3 + 0]);
	minYpLoc[1] = real(m->xyz[left_cv*3 + 1]);
	minYpLoc[2] = real(m->xyz[left_cv*3 + 2]);
      }
    }
  }

  MPI_Datatype mpit;
  //get mpi datatype to send
  mpit = MPI_GetType(maxYp);

  //parallel sync
  MPI_Allreduce(MPI_IN_PLACE, &maxYp, 1, mpit, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &minYp, 1, mpit, MPI_MIN, MPI_COMM_WORLD);

  //TODO: we need to check for the process which holds the smallest/largest y+ values
  //and broadcast that to every process for this to be valid. For now, print everywhere
  std::cout << "Max Y+ Location: " << maxYpLoc[0] << " " << maxYpLoc[1] << " " << maxYpLoc[2] << std::endl;
  std::cout << "Min Y+ Location: " << minYpLoc[0] << " " << minYpLoc[1] << " " << minYpLoc[2] << std::endl;

  //set print formats
  fout.setf(std::ios::scientific);
  fout.precision(16);

  if(p->GetRank() == 0){
    fout << "Max Y+: " << maxYp << std::endl;
    fout << "Min Y+: " << minYp << std::endl;
  }

}

template <class Type>
void Forces<Type>::AllocateForcesMemory(Int num_bcs, Int num_bodies)
{
  this->num_bodies = num_bodies;
  this->num_bcs = num_bcs;
  this->bodies = new CompositeBody<Type>[num_bodies+1];
  Mesh<Type>* m = space->m;
  Int nbedge = m->GetNumBoundaryEdges();
  yp = cf = qdot = cp = NULL;
  yp = new Type[nbedge];
  cf = new Type[nbedge];
  qdot = new Type[nbedge];
  cp = new Type[nbedge];
  //cumulate surface areas of each bc... summed abs value
  surfArea = new Type[(num_bcs+1)*3];

  if(surfArea == NULL || cp == NULL || qdot == NULL || cf == NULL || yp == NULL){
    std::cerr << "Solution memory allocation error -- Forces::AllocateForcesMemory()!!" 
	       << std::endl;
  }
  forcesAlloc = 1;
}

template <class Type>
void Forces<Type>::KillForcesMemory()
{
  if(forcesAlloc == 1){
    delete [] surfArea;
    delete [] yp;
    delete [] cf;
    delete [] cp;
    delete [] qdot;    
    delete [] bodies; //bury the bodies...
  }
  forcesAlloc = false;
}
