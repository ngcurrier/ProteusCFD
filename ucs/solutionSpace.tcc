#include "bc.h"
#include "mesh.h"
#include "parallel.h"
#include "create_functions.h"
#include "jacobian.h"
#include "residual.h"
#include "solve.h"
#include "sensors.h"
#include "solutionField.h"
#include "forces.h"
#include "dataInfo.h"
#include "limiters.h"
#include "vars_id.h"
#include "eqnset.h"
#include "param.h"
#include "turb.h"
#include "exceptions.h"
#include "customics.h"
#include "gaussian.h"
#include "postForce.h"
#include <sstream>
#include <fstream>
#include <iostream>

template <class Type>
SolutionSpace<Type>::SolutionSpace(Param<Type>* param, PObj<Real>* p, std::string name, 
				   TemporalControl<Real>& temporalControl):
  SolutionSpaceBase<Type>(name, temporalControl), param(param), p(p), turb(NULL), gaussian(NULL),
  q(NULL), qold(NULL), qoldm1(NULL), qgrad(NULL)
  
{
  isCopy = false;
  Int ierr = 0;
  this->timers.InitList(9);
  this->timers.CreateTimer("MapsBuildTimer");
  this->timers.CreateTimer("WallDistanceTimer");
  this->timers.CreateTimer("IterationTimer");
  this->timers.CreateTimer("ResidualTimer");
  this->timers.CreateTimer("LinearSolveTimer");
  this->timers.CreateTimer("JacobianAssembleTimer");
  this->timers.CreateTimer("GradientTimer");
  this->timers.CreateTimer("SolutionUpdateTimer");
  this->timers.CreateTimer("TurbulenceModelTimer");

  m = new Mesh<Type>;
  bc = new BoundaryConditions<Type>;
  m->SetParallelPointer(p);
  
  //assume that all mesh names are given w/o the full path, build the path here
  ierr = m->ReadPartedMesh(param->path+param->spacename);
  if(param->scaleMesh){
    m->ScaleMesh(1.0/param->ref_length);
  }
  if(ierr){
    Abort << "Partitioned mesh read failed";
  }

  //build the mesh maps.. this creates all ghost/phantom nodes as appropriate
  //needs to be done before memInit which happens after eqnset object creation
  this->timers.StartTimer("MapsBuildTimer");
  if(param->reorder == 1 && m->IsReordered() == false){
    //this builds only the maps up to psp()
    ierr = m->BuildMapsDecomp();
    //use reverse CuthillMcKee
    ierr = m->ReorderMeshCuthillMcKee(1);
    m->CleanMapsDecomp();
    ierr = m->ReorderC2nMap();
    //this builds all the maps with the new ordering
    ierr = p->ReorderCommNodes(m);
    ierr = m->BuildMaps();
    //the reordering messes up the partition file relations
    //just re-write them to fix the issue
    m->WriteParallelMesh(param->path+param->spacename);
  }
  else{
    ierr = m->BuildMaps();
  }
  if(ierr){
    Abort << "Mesh maps building failed";;
  }
  //clear old solution data from file if not restarting
  if(!param->useRestart){
    ClearSolutionFromFile();
  }

  //after reading bc file all internal counters etc. are set
  //the only exception is the qref states not read from file
  //explicitly do not yet point at eqnset->qref
  ierr = bc->ReadFile(param->path+param->spacename, m->GetMaximumFactag());
  if(ierr){
    Abort << "BC file read failed";
  }
  this->timers.StopTimer("MapsBuildTimer");

  //create eqnset object of correct type
  ierr = CreateEqnSet(this);
  if(ierr){
    Abort << "EqnSet creation failed";
  }

  //initialize bcobjs pointers, etc. which were not read from file
  //also sets the internal pointer to qref in the eqnset
  ierr = bc->InitInternals(eqnset);
  if(ierr){
    Abort << "EqnSet init internals failed";
  }

  //this finishes all of the internal initialization not related to the mesh (which
  //might be copied later
  Init();

  //Write solution of zeroth time step to check BCs, initialization, etc.
  WriteSolution();
}

template <class Type>
SolutionSpace<Type>::~SolutionSpace()
{
  for(UInt i = 0; i < fields.size(); i++){
    delete fields[i];
  }

  if(!isCopy){
    delete bc;
  }
  else{
    delete p;
  }

  delete param;
  delete m;
  delete eqnset;
  delete sensors;
  delete turb;
  delete forces;
  delete limiter;
  delete crs;
  delete grad;

  delete ceqnset;
  delete cparam;
  delete gaussian;

  CloseOutFiles();

  return;
}

template <class Type>
void SolutionSpace<Type>::Init()
{
  Int ierr = 0;

  //we only want open output files on the original solution space
  OpenOutFiles();

  //finish up the initialization of the PObj
  ierr = p->BuildCommMaps(m);
  if(ierr){
    Abort << "Mesh communication maps building failed";
  }
  ierr = m->CalcMetrics();
  if(ierr){
    Abort << "Mesh metrics calulation failed";
  }

  //init eqnset for solve - init solution vector storage (q, qold, qoldm1, qgrad, etc)
  eqnset->InitEqnSet();

  //set status flags for boundary conditions
  SetBCStatFlags(m, bc);

  //this logic is specific to error transport equations for error indicator
  if(param->viscous && param->errorTransport){
    DataInfo ETEViscData(eqnset->neqn, std::string("ETEVisc"));
    const std::vector<std::string>& names = eqnset->idata->GetNames();
    for(Int i = 0; i < names.size(); i++){
      std::string dofName = names[i];
      if(eqnset->idata->DofIsScalar(i)){
	ETEViscData.AddScalar(i, dofName+"_ViscErrorSrc");
      }
      else{
	ETEViscData.AddVector(i, dofName+"_ViscErrorSrc");
	i += 2;
      }
    }
    ETEViscData.Verify();
    AddField(ETEViscData, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  }

  //add field for gaussian source
  if(param->gaussianSource){
    AddField("gaussian");
  }

  std::cout << "INIT: gaussian" << std::endl;
  if(param->gaussianSource){
    gaussian = new GaussianSource<Type>(param->gaussianEqn, param->gaussianBCid, param->gaussianXloc, 
					param->gaussianYloc, param->gaussianAmpl, *this);
  }

  std::cout << "INIT: gradient" << std::endl;
  //allocate gradient object
  std::vector<Int> gradientList;
  Int nterms = eqnset->GetGradientsLocation(gradientList);
  Int* list = new Int[nterms];
  for(Int i = 0; i < nterms; i++){
    list[i] = gradientList[i];
  }
  Int weighted = true;
  grad = new Gradient<Type>(nterms, eqnset->neqn+eqnset->nauxvars, list, q, this, 
			    param->gradType, qgrad, weighted);
  delete [] list;

  std::cout << "INIT: LSQ coefficient calculation" << std::endl;
  //compute gradient coefficients once
  ComputeNodeLSQCoefficients(this);

  std::cout << "INIT: Limiter" << std::endl;
  //allocate memory required for our limiters, this lives here because it should be with the data
  limiter = new Limiter<Type>(this, q, qgrad, eqnset->neqn, eqnset->neqn+eqnset->nauxvars, grad->GetNterms(), "variableQ");

  std::cout << "INIT: sensors" << std::endl;
  //allocate sensors object if appropriate
  Int nnode = m->GetNumNodes();
  sensors = new Sensors<Type>(param->path+param->spacename+".sensors", eqnset->neqn, 
			      eqnset->nauxvars+eqnset->neqn, q, 
			      p, m->ipsp, m->psp, m->xyz, nnode);
  //search for sensor points once
  sensors->Search();
  
  std::cout << "INIT: forces" << std::endl;
  //allocate forces object if appropriate
  forces = new Forces<Type>(this);
  //allocate memory for composite surfaces
  forces->AllocateForcesMemory(bc->largest_bc_id, bc->num_bodies);
  //read in bc file again to look for the composite bodies
  BodiesReadFile(forces->bodies, param->path+param->spacename);

  //add field for dt local
  AddField("timestep");
  //add field for preconditioning term beta
  AddField("beta");
  SolutionField<Type>& betaF = GetField("beta");
  //set preconditioner beta equal to Ma^2 if below sonic
  Type V = param->GetVelocity(0);
  Type Mach = V;
  if(real(Mach) < 1.0){
    Type beta = MAX(real(param->betaMin), real(Mach*Mach));
    betaF.Fill(beta);
  }
  else{
    //turn off preconditioning
    betaF.Fill(1.0);
  }
  //add field for partition id
  AddField("partition");
  Type* part = GetFieldData("partition", FIELDS::STATE_NONE);
  Int rank = p->GetRank();
  for(Int i = 0; i < nnode; i++){
    part[i] = rank; 
  }
  p->UpdateGeneralVectors(part, 1);

  //add field for wall distance
  if(param->viscous){
    AddField("wallDistance");
    //allocate turbulence model if appropriate, has to be registered first
    //or read restart won't work
    CreateTurbModel(&turb, this);
  }


  if(param->useRestart){
    ReadRestartFile();
  }
  else{
    eqnset->SetInitialConditions();
  }
  //set custom ICs if we are not restarting and they are asked for
  //call this last in case we need to modify the conditions given by default
  if(param->customIcId && !param->useRestart){
    InitIC(param->customIcId, this);
  }
  //set all old q's to steady state if not restarting
  if(param->torder && !param->useRestart){
    Int neqn = eqnset->neqn;
    Int nvars = neqn + eqnset->nauxvars;
    memcpy(qold, q, sizeof(Type)*nnode*nvars);
    memcpy(qoldm1, q, sizeof(Type)*nnode*nvars);
    for(Int i = 0; i < nnode; i++){
      eqnset->NativeToConservative(&qold[i*nvars]);
      eqnset->NativeToConservative(&qoldm1[i*nvars]);
    }
  }
  p->UpdateGeneralVectors(q, eqnset->neqn+eqnset->nauxvars);

  //sanity check for parallel comms based on passing node coordinates
  ierr = p->CheckSanityCoords(m->GetNodeCoords());
  if(ierr){
    Abort << "Parallel sanity check of coordinates failed";
  }

  //if we require a complex eqnset and chem model for jacobian use
  //create them now, we do not allow copies of solution spaces to
  //own complex eqnsets b/c in most cases this will fail anyway
  //currently the only copies we use are complex anyway, CTSE does not
  //work in those situations
  if((param->requireComplex || param->mode != 0) && !isCopy){
    std::cout << "Param file requires creation of complex equation set OR mode is NOT zero, creating" << std::endl;
    cparam = new Param<RCmplx>(*param);
    ierr = CreateEqnSet(&ceqnset, cparam);
    if(ierr){
      Abort << "Creation of complex EqnSet object failed";
    }
  }
  else{
    cparam = NULL;
    ceqnset = NULL;
  }

  if(param->viscous){
    std::cout << "Wall distance results: " << std::endl;
    std::cout << "=======================" << std::endl;
    this->timers.StartTimer("WallDistanceTimer");
    Type* dist = GetFieldData("wallDistance", FIELDS::STATE_NONE);
    if(false){
      ComputeWallDist(dist, this);
    }
    else{
      ComputeWallDistOct(dist, this);
    }
    //Update dist vectors via MPI
    p->UpdateGeneralVectors(dist, 1);
    this->timers.StopTimer("WallDistanceTimer");
    this->timers.PrintTimer("WallDistanceTimer", timerOutFile);
    std::cout << std::endl << std::endl;
  }

  //this has to come after wall distance computation in some cases
  if(param->viscous){
    //this is restart respecting
    turb->Initialize();
  }

  //allocate field for adjoint variables if running in that mode
  if(param->mode == Adjoint){
    DataInfo AdjointData(eqnset->neqn, std::string("Adjoint"));
    const std::vector<std::string>& names = eqnset->idata->GetNames();
    for(Int i = 0; i < eqnset->neqn; i++){
      std::string dofName = names[i];
      if(eqnset->idata->DofIsScalar(i)){
	AdjointData.AddScalar(i, dofName+"_Adjoint");
      }
      else{
	AdjointData.AddVector(i, dofName+"_Adjoint");
	i += 2;
      }
    }
    AdjointData.Verify();
    AddField(AdjointData, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
  }

  //allocate jacobian memory for our Mesh object
  //this comes after walldistance search b/c of memory usage issues
  InitCRSSystem();
  //print out some info on CRS data size if implicit
  if(crs != NULL){
    crs->GetSystemMemUsage(true);
  }

  //pre-run gaussian to setup src distribution, needed for some BCs
  if(param->gaussianSource){
    gaussian->ApplyToResidual();
  }

  //add in the plugins
  postPlugins.push_back(new PostForcePlugin<Type>(*this));
}

//add a field by name alone - single scalar
template <class Type>
void SolutionSpace<Type>::AddField(std::string name)
{
  DataInfo idata(1, name);
  idata.AddScalar(0, name);
  idata.Verify();

  AddField(idata, FIELDS::STATE_NONE, FIELDS::VAR_INTERIOR);
}

//add a field with state information, a data descriptor, and location info
template <class Type>
void SolutionSpace<Type>::AddField(DataInfo dataInfo, Int stateType, Int varLocation)
{
  std::cout << "SolutionSpace: adding field " << dataInfo.GetName() << std::endl;
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type>& field = **it;
    if(field.IsNamed(dataInfo.GetName())){
      Abort << "SolutionSpace::AddField() Field of name " + dataInfo.GetName()
	+ " already exists... must have unique name";
      return;
    }
  }
  SolutionField<Type> * temp =  new SolutionField<Type>(*m, dataInfo,  stateType, varLocation);
  fields.push_back(temp);
}

//remove a field by name alone
template <class Type>
void SolutionSpace<Type>::RemoveField(std::string name)
{
  std::cout << "SolutionSpace: removing field " << name << std::endl;
  for(Int i = 0; i < fields.size(); i++){
    if(fields[i]->IsNamed(name)){
      delete fields[i];
      fields.erase(fields.begin()+i);
    }
  }
  Abort << "Field of name " + name + " does not exist... cannot remove";
}

//get a field by name alone
template <class Type>
SolutionField<Type> & SolutionSpace<Type>::GetField(std::string name)
{
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type> & field = **it;
    if(field.IsNamed(name)){
      return (field);
    }
  }
  Abort << "Solution field \'" + name + "\' not found! ";
  return (**fields.begin());
}

//get a const field by name alone
template <class Type>
const SolutionField<Type> & SolutionSpace<Type>::GetField(std::string name) const
{
  for(typename std::vector<SolutionField<Type>*>::const_iterator it = fields.begin(); it != fields.end(); ++it){
    const SolutionField<Type> & field = **it;
    if(field.IsNamed(name)){
      return (field);
    }
  }
  Abort << "Solution field \'" + name + "\' not found! ";
  return (**fields.begin()); //this will never be executed - squash compiler warnings
}


//get a field array by name and temporal state, returns a raw pointer
template <class Type>
Type* SolutionSpace<Type>::GetFieldData(std::string name, Int state)
{
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type> & field = **it;
    if(field.IsNamed(name)){
      return (field.GetData(state));
    }
  }
  Abort << "Solution field of name " + name + " does not exist... returning NULL";
  return NULL;
}

//check if a field exists by name
template <class Type>
Bool SolutionSpace<Type>::CheckField(std::string name) const
{
  for(typename std::vector<SolutionField<Type>*>::const_iterator it = fields.begin(); it != fields.end(); ++it){
    const SolutionField<Type> & field = **it;
    if(field.IsNamed(name)){
      return true;
    }
  }
  return false;
}

//writes a list of all available fields to output
template <class Type>
void SolutionSpace<Type>::WriteAvailableFields() const
{
  std::cout << "Registered fields on solution space " << this->name << ": " << std::endl;
  std::cout << "============================================================" << std::endl;
  for(typename std::vector<SolutionField<Type>*>::const_iterator it = fields.begin(); it != fields.end(); ++it){
    const SolutionField<Type>& field = **it;
    std::cout << "\t+ " << field.GetName() << std::endl;
  }
}

//function to validate fields requested from input parameters to be valid names on solution space
template <class Type>
void SolutionSpace<Type>::ValidateRequestedFields() const
{
  std::stringstream ss;
  Bool err = false;
  for(typename std::vector<std::string>::iterator it = param->fieldsRequested.begin(); 
      it != param->fieldsRequested.end(); ++it){
    std::string & fieldname = *it;
    if(!CheckField(fieldname)){
      ss << "Field requested: " << fieldname << " does not exist\n";
      err = true;
    }
  }
  if(err){
    Abort << ss.str();
  }
}

//performs all solution space operations required prior to solution, called only once
template <class Type>
void SolutionSpace<Type>::PreIterate()
{
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();

  //verify that all fields that are requested are actually present
  ValidateRequestedFields();

  //compute surface areas for BCObjs 
  //this should be redone anytime we alter the surface mesh...
  //maybe here is not the best place but for now??
  ComputeSurfaceAreas(this, 1);

  //if we are restarting, allow the restart file to set the iteration count
  if(param->useRestart){
    if(!param->preserveRestartCounters){
      this->iter = 1;
    }
  }
  else{ //not using restart
    // we write out the time information for xyz/vol if restarting
    // no need to set it, in fact, it destroys time accuracy if we do
    //we need to save the xyz locations of the nodes at t = 0
    memcpy(m->xyz_base, m->xyz, sizeof(Type)*(nnode+gnode)*3);
    //we also need to zero the displacements over time since we use
    //higher order accurate grid speed terms
    memcpy(m->xyzold, m->xyz, sizeof(Type)*(nnode*3));
    memcpy(m->xyzoldm1, m->xyz, sizeof(Type)*(nnode*3));
    //we need to set our gcl volumes to be equal to the current time
    //this is taken care of in restart routines if using restarted solution
    memcpy(m->vololdm1, m->vol, sizeof(Type)*nnode);
    memcpy(m->volold, m->vol, sizeof(Type)*nnode);
  }

  //this has to be precomputed so temporal jacobians
  //have timestep information
  Type residualDelta = 9999.0;
  if(real(residualnm1) > 1.0e-21){ // protect against FPE
    residualDelta = (residualnm1 - residual)/residualnm1;
  }
  param->UpdateCFL(this->iter, real(residualDelta));
  dtmin = ComputeTimesteps(this);

  this->timers.StartTimer("IterationTimer");

  residual = residualnm1 = 0.0;

  //this is here in case we have a strided jacobian update period and we are
  //restarting, otherwise the jacobians are not computed on the first step
  if((param->nSgs > 0) && param->useRestart){
    if(p->GetRank() == 0){
      this->timers.StartAccumulate("JacobianAssembleTimer");
    }
    ComputeJacobians(this);
    if(p->GetRank() == 0){
      this->timers.PauseAccumulateAndPrint("JacobianAssembleTimer", timerOutFile);
    }
    Abort.CheckForSoftAbort("Soft abort post jacobians");
  }
}

//performs operations required prior to Newton iterations, once per timestep
template <class Type>
void SolutionSpace<Type>::PreTimeAdvance()
{
  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;
  Int nnode = m->GetNumNodes();

  std::cout << "\n" << this->iter << ": --------------------------------------------------------------------------\n\n";
  std::cerr << "\n" << this->iter << ": --------------------------------------------------------------------------\n\n";

  if(param->torder){
    for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
      SolutionField<Type> & field = **it;
      field.EvolveInTime();
    }
    //we only use these q's for unsteady simulations so we need to ensure
    //that these can be contributed directly to the RHS of the unsteady residual
    for(Int i = 0; i < nnode; i++){
      eqnset->NativeToConservative(&qold[i*nvars]);
    }
  }
  if((param->nSgs > 0) && (this->iter % param->jacobianFreq == 0 || this->iter == 1)){
    if(p->GetRank() == 0){
      this->timers.StartAccumulate("JacobianAssembleTimer");
    }
    std::cout << "Computing jacobians" << std::endl;
    std::cerr << "Computing jacobians" << std::endl;
    ComputeJacobians(this);
    if(p->GetRank() == 0){
      this->timers.PauseAccumulateAndPrint("JacobianAssembleTimer", timerOutFile);
    }
    Abort.CheckForSoftAbort("Soft abort post jacobians");
  }
  //we need to save the movement data for the h.o. temporal accuracy here
  if(param->movement){
    //copy down the old displacements so we get h.o. accurate grid speeds
    memcpy(m->xyzoldm1, m->xyzold, sizeof(Type)*3*nnode);
    memcpy(m->xyzold, m->xyz, sizeof(Type)*3*nnode);
    //if we are using a modifiable mesh we need to update the GCL contributions
    //for this we need the volumes from two previous steps
    memcpy(m->vololdm1, m->volold, sizeof(Type)*nnode);
    memcpy(m->volold, m->vol, sizeof(Type)*nnode);
    if(param->viscous){
      std::cout << "Wall distance results: " << std::endl;
      std::cout << "=======================" << std::endl;
      this->timers.StartTimer("WallDistanceTimer"); 
      Type* dist = GetFieldData("wallDistance", FIELDS::STATE_NONE);
      //ComputeWallDist(dist, this);
      ComputeWallDistOct(dist, this);
      //Update dist vectors via MPI
      p->UpdateGeneralVectors(dist, 1);
      this->timers.StopTimer("WallDistanceTimer");
      std::cout << std::endl << std::endl;
    }
  }

  //calculate timesteps, do once per timestep
  Type residualDelta = 9999.0;
  if(real(residualnm1) > 1.0e-21){ // protect against FPE
    residualDelta = (residualnm1 - residual)/residualnm1;
  }
  param->UpdateCFL(this->iter, real(residualDelta));
  dtmin = ComputeTimesteps(this);
  
}

//performs solution iteration, called multiple times per timestep
template <class Type>
bool SolutionSpace<Type>::NewtonIterate()
{
  Int nnode = m->GetNumNodes();
  Int nbnode = m->GetNumBoundaryNodes();
  Int gnode = m->GetNumParallelNodes();

  Int neqn = eqnset->neqn;
  Int nvars = neqn + eqnset->nauxvars;

  Type dqNorm, dqNormGlobal;
  Int resnode;
  Type tempmax, tempold, deltaDq;
 
  if(param->movement){
    std::cout << "Moving mesh" << std::endl;
    std::cerr << "Moving mesh" << std::endl;
    MoveMesh(this);
  }

  //set BCs for next iteration
  std::cout << "Updating boundary conditions" << std::endl;
  std::cerr << "Updating boundary conditions" << std::endl;
  UpdateBCs(this);

  //Update solution vectors via MPI
  p->UpdateGeneralVectors(q, nvars);

  if(p->GetRank() == 0){
    this->timers.StartAccumulate("GradientTimer");
  }

  //if doing higher order or viscous compute gradients
  if((param->sorder > 1 && (this->iter > param->nFirstOrderSteps)) || param->viscous){
    std::cout << "Computing gradient" << std::endl;
    std::cerr << "Computing gradient" << std::endl;
    grad->Compute();
  }

  //run limiter computation, we need to limit the extrapolation so we have to convert
  //the Q vector to the appropriate variable set first
  if((this->iter%param->limiterRefresh == 0) && (param->sorder > 1) && param->limiter
     && (this->iter > param->nFirstOrderSteps) && (this->iter < param->limiterFreeze)){
    std::cout << "Computing limiters" << std::endl;
    std::cerr << "Computing limiters" << std::endl;
    limiter->Compute(this);
  }
  
  if(p->GetRank() == 0){
    this->timers.PauseAccumulateAndPrint("GradientTimer", timerOutFile);
  }

  //calculate residuals
  if(p->GetRank() == 0){
    this->timers.StartAccumulate("ResidualTimer");
  }
  std::cout << "Computing fluxes/residuals" << std::endl;
  std::cerr << "Computing fluxes/residuals" << std::endl;
  std::vector<Type> residualComponents = ComputeResiduals(this);
  residGlobal = residualComponents[0]; //first value is the total residual
  Abort.CheckForSoftAbort("Soft abort post residuals");
  if(p->GetRank() == 0){
    this->timers.PauseAccumulateAndPrint("ResidualTimer", timerOutFile);
  }
  
  //copy residual for storage
  residualnm1 = residual;
  residual = residGlobal;
  
  std::cout << "dt_min: " << dtmin;
  std::cout << " time: " << this->time;
  
  //find global minimum timestep
  Real dtmin_global = real(dtmin);
  MPI_Datatype mpit = MPI_GetType(dtmin_global);
  MPI_Allreduce(MPI_IN_PLACE, &dtmin_global, 1, mpit, MPI_MIN, MPI_COMM_WORLD);

  if(p->GetRank() == 0){
    residOutFile << std::endl << this->iter << "\t" <<  residGlobal << "\t";
    for(Int i = 0; i < neqn; i++){
      residOutFile << residualComponents[i+1] << "\t";
    }
  }

  if(param->nSgs > 0){
    std::cout << "\nComputing implicit solution" << std::endl;
    std::cerr << "\nComputing implicit solution" << std::endl;
    if(p->GetRank() == 0){
      this->timers.StartAccumulate("LinearSolveTimer");
    }
    crs->A->PrepareSGS();
    //set initial guess to zero
    crs->BlankX();
    deltaDq = crs->SGS(param->nSgs, NULL, NULL, NULL);
    //deltaDq = crs->GMRES(1, param->nSgs, 2, NULL, NULL, NULL);
#if 0
    CRSMatrix<Type> Store;
    Store.CopyMatrixStructure(crs->A);
    Type* xstore;
    Type* bstore;
    crs->AllocateVector(&xstore);
    crs->AllocateVector(&bstore);
    memcpy(bstore, crs->b, sizeof(Type)*(neqn*(nnode+gnode)));
    
    deltaDq = crs->GMRES(1, 10, 1, &Store, xstore, bstore);
    std::cout << "\nD||diag||: " << deltaDq;
    deltaDq = crs->GMRES(1, 10, 2, &Store, xstore, bstore);
    std::cout << "\nD||block-diag||: " << deltaDq;
    deltaDq = crs->GMRES(1, 10, 3, &Store, xstore, bstore);
    std::cout << "\nD||ILU0||: " << deltaDq;
    deltaDq = crs->GMRES(1, 10, 4, &Store, xstore, bstore);
    std::cout << "\nD||SGS||: " << deltaDq;
    crs->A->PrepareSGS();
    deltaDq = crs->SGS(10, NULL, NULL, NULL, true);
    std::cout << "\nD||sgs-dq||: " << deltaDq;
#endif
    std::cout << "D||sgs-dq||: " << deltaDq;
    if(p->GetRank() == 0){
      this->timers.PauseAccumulateAndPrint("LinearSolveTimer", timerOutFile);
    }  
  }
  else{
    std::cout << "Computing explicit solution" << std::endl;
    std::cerr << "Computing explicit solution" << std::endl;
    ExplicitSolve(*this);
  }
  Abort.CheckForSoftAbort("Soft abort post linear solve");

  if(p->GetRank() == 0){
    this->timers.StartAccumulate("SolutionUpdateTimer");
  }

  //check for infinite updates and if found zero the node's update
  for(Int j = 0; j < nnode; j++){
    bool zero = false;
    for(Int k = 0; k < neqn; k++){
      if(std::isnan(real(crs->x[j*neqn +k])) || std::isinf(real(crs->x[j*neqn+k]))){
	crs->x[j*neqn + k] = 0.0;
	std::stringstream ss;
	ss << "Update dQ[" << k << "] isnan OR isinf at " <<
	  m->xyz[j*3 + 0] << "  " << m->xyz[j*3 + 1] << " " << m->xyz[j*3 + 2] <<
	  " Zeroing dQ update for node\n";
	std::cerr << ss.str();
	zero = true;
      }
    }
    if(zero){
      for(Int k = 0; k < neqn; k++){
	crs->x[j*neqn + k] = 0.0;
      }
    }
  }

  for(Int j = 0; j < nnode; j++){
    eqnset->ApplyDQ(&crs->x[j*neqn], &q[j*nvars], &m->xyz[j*3]);
  }  

  for(Int j = 0; j < nnode; j++){
    for(Int k = 0; k < neqn; k++){
      if(std::isnan(real(q[j*nvars + k])) || std::isinf(real(q[j*nvars+k]))){
	std::stringstream ss;
	ss << "Updated Q[" << k << "] isnan OR isinf at " <<
	  m->xyz[j*3 + 0] << "  " << m->xyz[j*3 + 1] << " " << m->xyz[j*3 + 2];
	Abort << ss.str();
      }
    }
  }

  dqNorm = VecL2Norm(crs->x, (nnode+gnode)*neqn);
  dqNormGlobal = ParallelL2Norm(p, crs->x, nnode*neqn);

  if(p->GetRank() == 0){
    this->timers.PauseAccumulateAndPrint("SolutionUpdateTimer", timerOutFile);
  }  
  
  std::cout << "\t||dq||: " << dqNorm;
  
  if(p->GetRank() == 0){
    residOutFile << dqNormGlobal << "\t";
  }
  
  //find max residual and cv id
  tempmax = VecL2Norm(&crs->b[0], neqn);
  resnode = 0;
  tempold = tempmax;
  for(Int j = 1; j < nnode; j++){
    tempmax = MAX(tempmax, VecL2Norm(&crs->b[j*neqn], neqn));
    if(tempold != tempmax){
      tempold = tempmax;
      resnode = j;
    }
  }
  std::cout << "\tmax||R||: " << tempmax << " @node " << resnode;
  
  //find max dq and cv id
  tempmax = VecL2Norm(&crs->x[0], neqn);
  Int qnode = 0;
  tempold = tempmax;
  for(Int j = 1; j < nnode; j++){
    tempmax = MAX(tempmax, VecL2Norm(&crs->x[j*neqn], neqn));
    if(tempold != tempmax){
      tempold = tempmax;
      qnode = j;
    }
  }
  std::cout << "\tmax||dq||: " << tempmax << " @node " << qnode << std::endl;     
  
  //Update solution vectors via MPI
  p->UpdateGeneralVectors(q, nvars);

  if(p->GetRank() == 0){
    this->timers.StartAccumulate("TurbulenceModelTimer");
  }
  //update turbulence model
  if(param->viscous){
    turb->Compute();
    if(param->gcl){
      std::cout << "COMPUTE GCL FOR TURB MODEL!" << std::endl;
    }
  }
  if(p->GetRank() == 0){
    this->timers.PauseAccumulateAndPrint("TurbulenceModelTimer", timerOutFile);
  }  

  Type currentCFL = param->GetCFL();
  std::cout << "CFL: " << currentCFL << "\n";
  if(p->GetRank() == 0){
    residOutFile << currentCFL << "\t" ;
    residOutFile << dtmin_global << "\t";
  }
  
  //we have to compute the forces inside the Newton loop
  //if doing fluid-structure interactions, we need cl, and cm
  forces->Compute();
  
  if(std::isnan(real(residGlobal)) || std::isinf(real(residGlobal))){
    std::stringstream ss;
    ss << "Solution residual divergent!! " << std::endl;
    ss << "Coordinates of max residual: " << m->xyz[resnode*3 + 0] 
       << " " << m->xyz[resnode*3 + 1] << " " << m->xyz[resnode*3 + 2] << std::endl;
    ss << "Coordinates of max dq: " << m->xyz[qnode*3 + 0] 
       << " " << m->xyz[qnode*3 + 1] << " " << m->xyz[qnode*3 + 2] << std::endl;
    Abort << ss.str();
    nanflag = 1;
  }

  if(real(residGlobal) < real(this->temporalControl.newtonConvergence)){
    return true;
  }
  else{
    return false;
  }

}

//performs operations required after Newton iterations, once per timestep
template <class Type>
void SolutionSpace<Type>::PostTimeAdvance()
{
  Int rank = p->GetRank();

  //increment the temporal counter since this is the new timestep, only do this once per time
  // not on every inner Newton loop
  this->time += param->dt;

  this->iter++;
  
  //sample sensor data
  sensors->Sample();
  //write out sensor data
  sensors->Write(param->path+param->spacename, this->iter);
  
  //Find forces and moments, cl, yp, cf, etc.
  forces->Compute();
  //report cf, yp, qdot, cl, forces, etc.
  forces->Report();

  //loop over all post processing plugins
  for(typename std::vector<PostPlugin<Type>*>::iterator it =
	postPlugins.begin(); it != postPlugins.end(); ++it){
    (*it)->Compute();
    (*it)->Report();
  }
  
  if(param->solutionWrite != 0){
    if(this->iter % param->solutionWrite == 0){
      WriteSolution();
    }
  }

  if(param->writeRestartStep != 0){
    if(this->iter % param->writeRestartStep == 0){
      WriteRestartFile();
    }
  }

  if(rank == 0){
    timerOutFile << this->iter << ":\n";
    this->timers.PrintSplit("IterationTimer", timerOutFile);
    //take care of the comm. timers
    p->timers.PrintAccumulate("CommTimer", timerOutFile);
    //reset the comm timer b/c this is used to track iteration comm count
    p->timers.ResetAccumulate("CommTimer");
  }
}

//returns a plugin pointer given a particular name
template <class Type>
PostPlugin<Type>* SolutionSpace<Type>::GetPlugin(std::string pluginName)
{
  //loop over all post processing plugins
  for(typename std::vector<PostPlugin<Type>*>::iterator it =
	postPlugins.begin(); it != postPlugins.end(); ++it){
    if ((*it)->isName(pluginName)) return (*it);
  }
  std::stringstream ss;
  ss << "Could not find plugin of name: " << pluginName << "\n";
  throw std::runtime_error(ss.str());
}


template <class Type>
void SolutionSpace<Type>::RefreshForParam()
{
  //these are calls which must be made so that post copy of a solution space
  //if a parameter is modified, the change gets passed along, for design
  eqnset->UpdateQinf();
  if(param->gaussianSource){
    gaussian->SetXY(param->gaussianXloc, param->gaussianYloc);
    gaussian->SetAmplitude(param->gaussianAmpl);
  }
}

template <class Type>
void SolutionSpace<Type>::OpenOutFiles()
{
  Int neqn = eqnset->neqn;

  //set residual print formats
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);
  std::string residualFileName = param->path+param->spacename + ".residual";

  if(isCopy){
    residualFileName = param->path+param->spacename + "_copy.residual";
  }

  if(m->p->GetRank() == 0){
    //if not restarting or restarting but wiping counters, open a new file
    if(!param->useRestart || !param->preserveRestartCounters){
      residOutFile.open(residualFileName.c_str());
    }
    //if restarting, append to the current file
    else{
      residOutFile.open(residualFileName.c_str(), std::ios::app);
    }
    if(!residOutFile.is_open()){
      std::stringstream ss;
      ss << "Residual file: " << residualFileName << " NOT opened!" << std::endl;
      Abort << ss.str();
    }
    //todo: print residual headers to file for plotting purposes
    residOutFile << "iter\tresGlobal\t";
    for(Int i = 0; i < neqn; i++){
      residOutFile << "resQ_" << i << "\t";
    }
    residOutFile << "||dq||\t" << "CFL\t" << "dtmin";
  }
  std::string timerFile = param->path+param->spacename + ".timer";
  if(isCopy){
    timerFile =  param->path+param->spacename + "_copy.timer";
  }
  if(m->p->GetRank() == 0){
    if(!param->useRestart || !param->preserveRestartCounters){
      timerOutFile.open(timerFile.c_str());
    }
    else{
      timerOutFile.open(timerFile.c_str(), std::ios::app);
    }
    if(!timerOutFile.is_open()){
      std::stringstream ss;
      ss << "Timer file: " << timerFile << " NOT opened!" << std::endl;
      Abort << ss.str();
    }
  }

  for(typename std::vector<PostPlugin<Type>*>::iterator it =
	postPlugins.begin(); it != postPlugins.end(); ++it){
    (*it)->WriteTabularHeader();
  }

}

template <class Type>
void SolutionSpace<Type>::CloseOutFiles()
{
  if(residOutFile.is_open()){
    residOutFile.close();
    timerOutFile.close();
  }
}

//initialize parallel compressed row storage for implicit solution
template <class Type>
void SolutionSpace<Type>::InitCRSSystem()
{
  if(eqnset == NULL){
    Abort << "Eqnset pointer is NULL in InitCRSSytem()... FAILING";
  }
  crs = new CRS<Type>;
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  crs->Init(nnode, gnode, eqnset->neqn, m->ipsp, m->psp, m->p);
}

template <class Type>
void SolutionSpace<Type>::WriteSolution()
{
  Int rank = p->GetRank();
  std::stringstream ss;
  ss.str("");
  ss << rank;
  std::string filename = param->path+param->spacename + "." + ss.str() + ".h5";
  //solution writing is turned off, likely running complex, etc.
  if(param->solutionWrite < 0){
    return;
  }
  std::cout << "-----------------------------------------------------------------" << std::endl;
  std::cout << "HDF_IO: Writing solution to file " << filename << std::endl;  
  std::cout << "-----------------------------------------------------------------" << std::endl;

  if(param->movement){
    m->WriteCurrentCoords(param->path+param->spacename, this->iter);
  }


  std::string directoryBase = "/Solution/";
  if(param->solutionTagStep){
    ss.str("");
    ss << this->iter;
    std::string timeflag = "timestep-" + ss.str() + "/";
    //modify path with timestep information if running unsteady
    directoryBase += timeflag;
    if(real(param->dt) > 0.0){
      hid_t h5out = HDF_OpenFile(filename, 1);
      if(h5out < 0){
	Abort << "SolutionSpace::WriteSolution() could not open file -- " + filename;
      }
      std::string timevaluedir = "/SolutionTime/" + timeflag + "/";
      HDF_WriteScalar(h5out, timevaluedir, "time", &this->time);
      H5Fclose(h5out);
    }
		    
  }

  //loop over the fields and write them to hdf5
  for(typename std::vector<std::string>::iterator it = param->fieldsRequested.begin(); 
      it != param->fieldsRequested.end(); ++it){
    std::string & fieldname = *it;
    SolutionField<Type>& field = GetField(fieldname);
    field.WriteH5(filename, directoryBase);
  }
}

template <class Type>
void SolutionSpace<Type>::ClearSolutionFromFile()
{
  std::stringstream ss;
  ss.str("");
  Int rank = p->GetRank();
  ss << rank;
  std::string filename = param->path+param->spacename + "." + ss.str() + ".h5";
  std::string rsfilename = param->path+param->spacename + "." + ss.str() + ".rs";

  //CUE MINI RANT:
  //  Why is HDF5 not clever enough to freaking delete space after a dataset
  //  no longer has any links to it, or even if I tell it... design FLAW!
  //  Okay, I feel better, just write the mesh back to a new file... 
  //  accomplishes what we want

  std::string rmline = "rm " + filename;
  system(rmline.c_str());
  rmline = "rm " + rsfilename;
  system(rmline.c_str());
  
  std::cout << "Clearing old solution data" << std::endl;
  std::cout << "Clearing old restart data" << std::endl;

  m->WriteParallelMesh(param->path+param->spacename);
}

template <class Type>
void SolutionSpace<Type>::WriteRestartFile()
{
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();
  std::ostringstream temposs;
  temposs.str("");
  temposs << p->GetRank();
  std::string filename = param->path+param->spacename + "." + temposs.str() + ".rs";
  
  std::cout << "RESTART WRITEFILE: writing restart file --> " << filename << std::endl;

  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type>& field = **it;
    if(field.IsTemporal()){
      field.WriteH5(filename, "/Fields/", FIELDS::STATE_TIME);
    }
  }

  hid_t h5out = HDF_OpenFile(filename, 1);
  HDF_WriteScalar(h5out, "/Solver State/", "Iteration Count", &this->iter);
  HDF_WriteScalar(h5out, "/Solver State/", "GCL State", &param->gcl);
  HDF_WriteScalar(h5out, "/Solver State/", "time", &this->time);
  HDF_WriteArray(h5out, "/Mesh State/", "Volume N+1", m->vol, nnode);
  HDF_WriteArray(h5out, "/Mesh State/", "Volume N", m->volold, nnode);
  HDF_WriteArray(h5out, "/Mesh State/", "Volume N-1", m->vololdm1, nnode);
  HDF_WriteArray(h5out, "/Mesh State/", "Nodal Coordinates N+1", m->xyz, 3*(nnode+gnode));
  HDF_WriteArray(h5out, "/Mesh State/", "Nodal Coordinates Base", m->xyz_base, 3*(nnode+gnode));
  HDF_WriteArray(h5out, "/Mesh State/", "Nodal Coordinates N", m->xyzold, 3*nnode);
  HDF_WriteArray(h5out, "/Mesh State/", "Nodal Coordinates N-1", m->xyzoldm1, 3*nnode);
  H5Fclose(h5out);
}

template <class Type>
void SolutionSpace<Type>::ReadRestartFile()
{
  std::ostringstream temposs;
  temposs.str("");
  temposs << p->GetRank();
  std::string filename = param->path+param->spacename + "." + temposs.str() + ".rs";
  Int nnode = m->GetNumNodes();
  Int gnode = m->GetNumParallelNodes();

  std::cout << "RESTART READFILE: reading restart file --> " << filename << std::endl;

  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type>& field = **it;
    if(field.IsTemporal()){
      field.ReadH5(filename, "/Fields/", FIELDS::STATE_TIME);
      std::cout << "RESTART READFILE: reading field " << field.GetName() << std::endl;
    }
  }

  hid_t h5in = HDF_OpenFile(filename, 0);
  HDF_ReadScalar(h5in, "/Solver State/", "Iteration Count", &this->iter);
  HDF_ReadScalar(h5in, "/Solver State/", "time", &this->time);
  Int gclcheck = 0;
  HDF_ReadScalar(h5in, "/Solver State/", "GCL State", &gclcheck);
  if(gclcheck != param->gcl){
    Abort << "WARNING: GCL state not sane across restart, this is not going to end well";
  }
  HDF_ReadArray(h5in, "/Mesh State/", "Volume N+1", &m->vol, &nnode);
  m->SetNumNodes(nnode);
  HDF_ReadArray(h5in, "/Mesh State/", "Volume N", &m->volold, &nnode);
  HDF_ReadArray(h5in, "/Mesh State/", "Volume N-1", &m->vololdm1, &nnode);
  Int ntnodes = 3*(nnode+gnode);
  HDF_ReadArray(h5in, "/Mesh State/", "Nodal Coordinates N+1", &m->xyz, &ntnodes);
  HDF_ReadArray(h5in, "/Mesh State/", "Nodal Coordinates Base", &m->xyz_base, &ntnodes);
  ntnodes = 3*nnode;
  HDF_ReadArray(h5in, "/Mesh State/", "Nodal Coordinates N", &m->xyzold, &ntnodes);
  HDF_ReadArray(h5in, "/Mesh State/", "Nodal Coordinates N-1", &m->xyzoldm1, &ntnodes);
  H5Fclose(h5in);
}


template <class Type>
void SolutionSpace<Type>::PrintTimers()
{
  if(p->GetRank() == 0){
    timerOutFile << "\nTimer Totals (Accumulated): " << std::endl;
    timerOutFile << "------------------------------" << std::endl;
    this->timers.PrintAllAccumulate(timerOutFile);
    this->timers.PrintAllTimers(timerOutFile);
    p->timers.PrintAccumulate("CommTotalTimer", timerOutFile);
  }
}


//copy constructor for solution space from one type (i.e. double) to another (i.e. complex)
//used mostly for sensitivity derivative calculations
template <class Type> template <class Type2>
SolutionSpace<Type>::SolutionSpace(const SolutionSpace<Type2>& spaceToCopy) : 
  SolutionSpaceBase<Type>(spaceToCopy.name, spaceToCopy.temporalControl), turb(NULL), gaussian(NULL)
{
  isCopy = true;

  //copy iteration counter since this controls second order switching and all sorts of other
  //internal goodies
  this->iter = spaceToCopy.iter;

  this->timers.InitList(9);
  this->timers.CreateTimer("MapsBuildTimer");
  this->timers.CreateTimer("WallDistanceTimer");
  this->timers.CreateTimer("IterationTimer");
  this->timers.CreateTimer("ResidualTimer");
  this->timers.CreateTimer("LinearSolveTimer");
  this->timers.CreateTimer("JacobianAssembleTimer");
  this->timers.CreateTimer("GradientTimer");
  this->timers.CreateTimer("SolutionUpdateTimer");
  this->timers.CreateTimer("TurbulenceModelTimer");

  Int ierr = 0;
  //this has to copy or create ALL of the relevant things in a solution space so that it can be used
  //as a replacement for the space passed in

  std::cout << "SOLUTION SPACE COPY: ParamCopy()" << std::endl;
  param = new Param<Type>(*spaceToCopy.param);
  p = new PObj<Type>();
  bc = spaceToCopy.bc;

  std::cout << "SOLUTION SPACE COPY: MeshCopy()" << std::endl;
  m = new Mesh<Type>(*spaceToCopy.m);
  m->SetParallelPointer(p);

  //create eqnset object of correct type
  std::cout << "SOLUTION SPACE COPY: CreateEqnSet()" << std::endl;
  ierr = CreateEqnSet(this);
  if(ierr){
    Abort << "EqnSet creation failed";
  }
  
  //this allocates memory, CRS structure, etc. internally on our copied solutionSpace
  std::cout << "SOLUTION SPACE COPY: calling Init()" << std::endl;
  Init();

  std::cout << "SOLUTION SPACE COPY: copy all fields" << std::endl;
  //now we need to duplicate the state of the solutionSpace by copying across all fields
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); 
      it != fields.end(); ++it){
    SolutionField<Type>& localField = **it;
    std::string name = localField.GetName();
    const SolutionField<Type2>& copiedField = spaceToCopy.GetField(name);
    localField = copiedField;
  }

  if(ierr){
    Abort << "SolutionSpace copy constructor fatal failure -- LATEZ!";
  }
}


