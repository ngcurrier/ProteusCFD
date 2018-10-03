#include "general.h"
#include "exceptions.h"
#include "bc.h"
#include "mesh.h"
#include "param.h"
#include "eqnset.h"
#include "create_functions.h"
#include "parallel.h"
#include "timer.h"
#include "customics.h"
#include "derivatives.h"
#include "portFileio.h"
#include "move.h"
#include "composite.h"
#include "solutionSpaceBase.h"
#include "solutionSpace.h"
#include "dataInfo.h"
#include "solve.h"
#include "solutionOperations.h"
#include "fluid_structure.h"
#include "temporalControl.h"
#include "pythonInterface.h"
#include <fenv.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <sys/utsname.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>

using namespace std;

int main(int argc, char* argv[]){

#ifdef _DEBUG
  //enable exceptions so we can trap NaNs, etc.
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

  Int rank, np;
  Int mode = 0;
  Int ndv = 0;
  Real* dObjdBeta = NULL;
  
  std::vector<Param<Real>* > paramList; 
  SolutionOrdering<Real> operations;
  TemporalControl<Real> temporalControl;

  string stdoutname = ".out.";
  string stderrname = ".err.";
  string tempname;
  stringstream temposs;

  TimerList timers(5);
  timers.CreateTimer("MPI_InitTimer");
  timers.CreateTimer("SolveTimer");
  timers.CreateTimer("DesignTimer");
  timers.CreateTimer("MovementTimer");

  //create parallel object for comms
  timers.StartTimer("MPI_InitTimer");
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  timers.StopTimer("MPI_InitTimer");

  PObj<Real> pobj;
  rank = pobj.GetRank();
  //remove Abort.Log if its present
  if(rank == 0){
    remove("Abort.Log");
  }
  
  if(!((argc == 2) || (argc == 3))){
    cerr << "Invalid arguments!!!" << endl;
    cerr << "USE: " << argv[0] << " <casename>" << endl;
    cerr << "OR" << endl;
    cerr << "USE: " << argv[0] << " <casename> <design type>" << endl;
    cerr << "<design type> - none = 0" << endl;
    cerr << "<design type> - objective f-n evaluation = 1" << endl;
    cerr << "<design type> - direct = 2" << endl;
    cerr << "<design type> - adjoint = 3" << endl;
    cerr << "<design type> - CTSE = 4" << endl;
    cerr << "<design type> - GRID SMOOTHING = 5" << endl;
    cerr << "<design type> - Compute Mesh Sensitivity = 6" << endl;
    cerr << "<design type> - Finite Difference = 7" << endl;

    Param<Real> tmp;
    tmp.PrintAllParams();

    MPI_Finalize();
    return (1);
  }

  std::string casestring = argv[1];
  size_t pos = casestring.rfind('/');
  std::string pathname;
  if(pos != std::string::npos){
    pathname = casestring.substr(0, pos+1);
    casestring = casestring.substr(pos);
  }
  else{
    pathname = "./";
  }

  //set pathname in abort class
  Abort.rootDirectory = pathname;

  if(argc == 3){
    temposs.clear();
    temposs.str("");
    temposs << argv[2];
    temposs >> mode;
  }

  //redirect std out
  temposs.clear();
  temposs.str("");
  temposs << rank;
  stdoutname = pathname + casestring + stdoutname + (temposs.str()); 
  stderrname = pathname + casestring + stderrname + (temposs.str());
  //only write to file if we are the last process
#ifndef _DEBUG
  if(pobj.GetRank() == pobj.GetNp()-1){
    freopen(stdoutname.c_str(), "w", stdout);
    freopen(stderrname.c_str(), "w", stderr);
  }
  else{
    freopen("/dev/null", "w", stdout);
    freopen("/dev/null", "w", stderr);
  }
#else
  //if _DEBUG is defined write all files
  freopen(stdoutname.c_str(), "w", stdout);
  freopen(stderrname.c_str(), "w", stderr);
#endif

  //get hostname, etc.
  struct utsname sysinfo;
  uname(&sysinfo);
  cout << "System Name: "<<sysinfo.sysname<<endl;
  cout << "Host Name: "<<sysinfo.nodename<<endl;
  cout << "Release(Kernel) Version: "<<sysinfo.release<<endl;
  cout << "Kernel Build Timestamp: "<<sysinfo.version<<endl;
  cout << "Machine Arch: "<<sysinfo.machine<<endl;
  cout << "Domain Name: "<<sysinfo.domainname<<endl;
  cout << "PID: " << (int)getpid() << endl;

  HDF_TurnOffErrorHandling();
  if(ReadParamFile(paramList, casestring, pathname)){
    Abort << "Error in param read";
    return (-1);
  }
  if(ReadSolutionOrdering(operations, casestring, pathname)){
    Abort << "Error in solution ordering read";
    return (-1);
  }
  if(ReadTemporalControl(temporalControl, casestring, pathname)){
    Abort << "Error in temporal control read";
    return (-1);
  }
  for(UInt i = 0; i < paramList.size(); ++i){
    paramList[i]->mode = mode;
  }
  //here if we are doing anything related to design, we must modify the parameters
  //according to the design file
  if(mode == ObjectiveEval || mode == GridSmoothing || mode == MeshSensitivity || 
     mode == Direct || mode == Adjoint || mode == CTSE || mode == FiniteDifference){
    for(UInt i = 0; i < paramList.size(); i++){
      paramList[i]->UpdateForDesign();
    }
  }
  //print out params after we have modified them as necessary
  for(std::vector<Param<Real>*>::iterator it = paramList.begin();
      it != paramList.end(); ++it){
    (*it)->PrintSolverParams();
  }
  //setup solution space container
  vector<SolutionSpaceBase<Real>*> solSpaces;
  for(std::vector<Param<Real>*>::iterator it = paramList.begin();
      it != paramList.end(); ++it){
    Param<Real>* param = *it;
    SolutionSpaceBase<Real>* solSpace;
    if(param->spacename == "structure"){
      solSpace = new STRUCTDYN::Structure(param, param->spacename, temporalControl);
    }
    else if(param->spacename == "typicalSection"){
      solSpace = new STRUCTDYN::TypicalSection(param, param->spacename, temporalControl);
    }
    else{
      solSpace = new SolutionSpace<Real>(param, &pobj, param->spacename, temporalControl);
    }
    solSpaces.push_back(solSpace);
  }
  for(std::vector<SolutionSpaceBase<Real>*>::iterator it = solSpaces.begin(); 
      it != solSpaces.end(); ++it){
    SolutionSpaceBase<Real> & space = **it;
    space.WriteAvailableFields();
  }
  //setup the solution operations (advanements, transfers, etc.)
  if(operations.Finalize(solSpaces)){
    MPI_Finalize();
    return(-1);
  }

  //call solver subsystem unless we are smoothing a mesh
  if(mode != GridSmoothing && mode != MeshSensitivity && mode != CTSE){
    timers.StartTimer("SolveTimer");
    Solve(solSpaces, operations);
    timers.StopTimer("SolveTimer");
  }

  //if mode is objective function evaluation
  if(mode == ObjectiveEval){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    Real ObjFunction = Compute_Obj_Function(space);
    string fullpath = space.param->path + space.param->spacename;
    cout << "\n\nOBJ-FUNCTION: " << ObjFunction << "\n\n" << endl;
    if(pobj.GetRank() == 0){
      WriteFunctionToDesignFile(fullpath, ObjFunction);
    }
  }

  //if mode is grid motion/smoothing
  if(mode == GridSmoothing){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    timers.StartTimer("MovementTimer");
    MoveMeshDesign(space.param->path + space.param->spacename, space.m, space.bc);
    timers.StopTimer("MovementTimer");

    //write the mesh  
    space.m->WriteParallelMesh(space.param->path+space.param->spacename);
  }

  //if mode is mesh sensitivity derivatives
  if(mode == MeshSensitivity){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    timers.StartTimer("MovementTimer");
    //compute the mesh sensitivities
    Compute_dXdB(space);
    //Compute_dXdB_FD(space);
    timers.StopTimer("MovementTimer");
  }
  
  std::cout.setf(std::ios::scientific);
  std::cout.precision(16);

  //if mode is any of the solution sensitivity functionality 
  if(mode == Direct || mode == Adjoint || mode == CTSE || mode == FiniteDifference){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    Param<Real>* param = space.param;
    timers.StartTimer("DesignTimer");
    //turn off restarting and solution writing at this point
    param->writeRestartStep = 0;
    param->useRestart = 0;
    param->solutionWrite = -1;
    
    //get number of design variables
    string fullpath = space.param->path + space.param->spacename;
    ndv = GetNdvDesignFile(fullpath);
    
    //NOTE: this completely displaces solution from above
    //must do this last!
    dObjdBeta = new Real[ndv];
    
    //Debugging output for checking dQdBetas
    //Compute_dQdBeta_FD(dQdB, eqnset, 0);
    //Real* dQdB_C = new Real[m.nnode*eqnset->neqn];
    //Compute_dQdBeta_CTSE(dQdB_C, eqnset, 0); 
    //for(Int i = 0; i < m.nnode; i++){
    //  for(Int j = 0; j < eqnset->neqn; j++){
    //    cout << dQdB_C[i*eqnset->neqn + j] << " " << dQdB[i*eqnset->neqn + j] << endl;
    //  }
    //}
  }

  //if mode is direct solution sensitivity
  if(mode == Direct){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    for(Int beta = 0; beta < ndv; beta++){
      dObjdBeta[beta] = Compute_dObjdBeta_Direct(space, beta);
      std::cout.setf(std::ios::scientific);
      std::cout.precision(16);
      cout << "dI/dBeta_direct[" << beta << "]: " << dObjdBeta[beta] << endl;
    } 
    for(Int beta = 0; beta < ndv; beta++){
      cout << "dI/dBeta_direct[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
  }

  //if mode is adjoint sensitivity
  if(mode == Adjoint){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    //adjoint
    Real* lambda = space.GetFieldData("Adjoint", FIELDS::STATE_NONE);
    //Compute_dObjdQ(lambda, eqnset, 0);
    Compute_Adjoint_Variables_II(lambda, space);
    for(Int beta = 0; beta < ndv; beta++){
      dObjdBeta[beta] = Compute_dObjdBeta_Adjoint(lambda, space, beta);
      std::cout.setf(std::ios::scientific);
      std::cout.precision(16);
      cout << "dI/dBeta_Adjoint[" << beta << "]: " << dObjdBeta[beta] << endl;
    }   
    for(Int beta = 0; beta < ndv; beta++){    
      cout << "dI/dBeta_Adjoint[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
  }

  //if mode is complex taylor series expansion sensitivity
  if(mode == CTSE){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    //complex differencing
    for(Int beta = 0; beta < ndv; beta++){
      dObjdBeta[beta] = Compute_dObjdBeta_CTSE(space, operations, beta);
      std::cout.setf(std::ios::scientific);
      std::cout.precision(16);
      cout << "dI/dBeta_CTSE[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
    for(Int beta = 0; beta < ndv; beta++){    
      cout << "dI/dBeta_CTSE[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
  }

  //if mode is finit difference sensitivity
  if(mode == FiniteDifference){
    Real h = 1.0e-8;
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    //finite differencing
    Real ObjFunctionBase = Compute_Obj_Function(space);

    Real* bounds = new Real[ndv*2];
    Int* dvType = new Int[ndv];
    Real* x = new Real[ndv];
    Real* grad = new Real[ndv];
    Real f;
    
    if(ReadDesignFile(space.param->path+space.param->spacename, &ndv, x, bounds, &f, grad, dvType)){
      Abort << "Could not read design file";
    }

    for(Int beta = 0; beta < ndv; beta++){
      Real xtemp = x[beta];
      //perturb beta, write to design file for movement
      x[beta] = h;
      if(pobj.GetRank() == 0){
	WriteDesignFile(space.param->path+space.param->spacename, ndv, x, bounds, f, grad, dvType);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      //move mesh then compute the objective function
      MoveMeshDesign(space.param->path + space.param->spacename, space.m, space.bc);
      Perturb_Param(beta, *space.param, 1);
      space.RefreshForParam();
      space.eqnset->SetInitialConditions();
      Solve(solSpaces, operations);
      Real ObjFunctionPerturbed = Compute_Obj_Function(space);
      dObjdBeta[beta] = (ObjFunctionPerturbed - ObjFunctionBase)/h;
      //reset the design variable for the next pass
      x[beta] = xtemp;
      //reset param file for next pass
      Perturb_Param(beta, *space.param, -1);
      space.RefreshForParam();
      //reset the grid for the next pass, keeps accumulation of movement from occurring
      Int nnode = space.m->GetNumNodes();
      Int gnode = space.m->GetNumParallelNodes();
      memcpy(space.m->xyz, space.m->xyz_base, sizeof(Real)*3*(nnode+gnode));
      space.m->CalcMetrics();
      std::cout.setf(std::ios::scientific);
      std::cout.precision(16);
      cout << "dI/dBeta_FD[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
    //make sure the design file is untouched as far as beta goes
    for(Int i = 0; i < ndv; i++){
      grad[i] = dObjdBeta[i];
    }
    if(pobj.GetRank() == 0){
      WriteDesignFile(space.param->path+space.param->spacename, ndv, x, bounds, f, grad, dvType);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(Int beta = 0; beta < ndv; beta++){    
      cout << "dI/dBeta_FD[" << beta << "]: " << dObjdBeta[beta] << endl;
    }
    delete [] bounds;
    delete [] dvType;
    delete [] x;
    delete [] grad;
  }

  //if any of the sensitivity modes, write output to file
  if(mode == Direct || mode == Adjoint || mode == CTSE || mode == FiniteDifference){
    SolutionSpace<Real>& space = *dynamic_cast<SolutionSpace<Real>*>(solSpaces[0]);
    timers.StopTimer("DesignTimer");
    string fullpath = space.param->path + space.param->spacename;
    if(pobj.GetRank() == 0){
      WriteGradToDesignFile(fullpath, dObjdBeta);
    }
    delete [] dObjdBeta;
  }
  
  //write solution at the end of the run always, unless we have
  //nothing to write out.
  if(mode != MeshSensitivity && mode != GridSmoothing){
    for(std::vector<SolutionSpaceBase<Real>*>::iterator it = solSpaces.begin(); 
	it != solSpaces.end(); ++it){
      SolutionSpaceBase<Real> & space = **it;
      space.WriteSolution();
      space.WriteRestartFile();
    }
  }

  //print all timers
  timers.PrintAllTimers(std::cout);

  //this is to clean up our solution spaces, etc.
  for(std::vector<SolutionSpaceBase<Real>*>::iterator it = solSpaces.begin(); it != solSpaces.end(); ++it){
    SolutionSpaceBase<Real>* space = *it;
    space->PrintTimers();
    delete space;
  }
  solSpaces.clear();

  MPI_Finalize();
  
  return (0);
}
