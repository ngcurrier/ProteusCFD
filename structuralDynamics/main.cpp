#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "timer.h"
#include "element_lib.h"
#include "structparam.h"
#include "utility.h"
#include "explicit.h"
#include "implicit.h"
#include "forces.h"
#include "bc.h"
#include "io.h"

using namespace std;
using namespace STRUCTDYN;

int main(int argc, char* argv[])
{
  int i, j, jj, t;
  stringstream ss;
  int mode;

  //set this flag to true if we want to hand
  //specify mass, stiffness, forcing, etc. in the config file
  //set this flag to false if using the Finite Element model
  bool specifiedMatrices = false;

  //time parameters
  double tfinal;
  //# of degrees of freedom
  int dof;
  //matrices (mass, damping, stiffness)
  double* m, *c, *k;
  //diagonalized mass matrix -- only for central and half-step central
  double* diagm = NULL;
  double* eldiagm = NULL;
  //newmark-beta parameters
  double gamma, beta;
  //wilson-theta parameters
  double theta;

  //time level information
  int nsteps;
  double dt;

  //solution memory
  double* x, *xd, *xdd;
  
  if(argc != 4){
    cerr << "USAGE: " << argv[0] << " t_final timestep mode" << endl;
    cerr << "modes:  0 - analytic (one degree of freedom only!)\n";
    cerr << "        1 - 2nd O RK - Raulston's method\n";
    cerr << "        2 - 2nd O RK - Trapezoidal rule\n";
    cerr << "        3 - 3rd O RK - Simpson's 1/3rd rule\n";
    cerr << "        4 - 4th O RK - Classical\n";
    cerr << "        5 - Classical central difference (full 2nd order - non self starting)\n";
    cerr << "        6 - Half step central difference (technically only 1st order) \n";
    cerr << "        7 - Newmark-Beta - average acceleration \n";
    cerr << "        8 - Newmark-Beta - linear acceleration \n";
    cerr << "        9 - Newmark-Beta - Fox-Goodwin         \n";
    cerr << "       10 - Wilson-Theta\n" << endl;
    return (-1);
  }
  ss.clear();
  ss.str("");
  ss << argv[1];
  ss >> tfinal;
  ss.clear();
  ss.str("");
  ss << argv[2];
  ss >> dt;
  ss.clear();
  ss.str("");
  ss << argv[3];
  ss >> mode;

  nsteps = (int)((double)tfinal/dt) + 1;

  //create forces object
  SForces forces;
  //create BCs object
  BC bcs;
  //parameters file, stores the mesh and case info
  SParam param;

  bcs.Read();
  bcs.SetParam(&param);

  if(specifiedMatrices){
    Read1DCase(&dof, nsteps, dt, &x, &xd, &xdd, &m, &c, &k, &forces);
  }
  else{
    std::string fileName = "infile.msh";
    dof = param.ReadMesh(fileName, &forces);
    bool dynamic = true;
    bool gravity = false;
    bool damping = true;
    bool diagmass;
    double alpha = 0.0002;
    double beta = 0.0003;

    forces.Init(param.nnodes, 0.0, param.dof);
    forces.Read();

    m = new double[dof*dof];
    c = new double[dof*dof];
    k = new double[dof*dof];

    //allocate solution memory
    x = new double[dof*(nsteps+1)];
    xd = new double[dof*(nsteps+1)];
    xdd = new double[dof*(nsteps+1)];

    //zero mass, stiffness and damping
    for(i = 0; i < dof*dof; i++){
      m[i] = k[i] = c[i] = 0.0;
    }
    
    param.SetICs(x, xd);

    int ndofpe = param.elements[0].GetElemDOF();
    double* elm = new double[ndofpe*ndofpe];
    double* els = new double[ndofpe*ndofpe];
    
    if(mode == 5 || mode == 6){
      //we can only use a diagonalized mass matrix for 
      //central difference and half-step central
      //will NOT work correctly for RK methods
      diagmass = true;
      diagm = new double[dof*dof];
      eldiagm = new double[ndofpe*ndofpe];
    }
    else{
      diagmass = false;
    }

    for(i = 0; i < param.nelem[0]; i++){
      param.elements[i].Compute(elm, els, dynamic, gravity, &param);
      if(diagmass){
	param.elements[i].LumpHRZ(eldiagm, &param);
      }
      param.elements[i].Transform(elm, eldiagm, els, dynamic, &param);
      param.elements[i].Assemble(elm, eldiagm, els, m, diagm, k, c, damping, 
				 gravity, alpha, beta, &param);
    }

    delete [] elm;
    delete [] eldiagm;
    delete [] els;
  }

  cout << "COMPUTED TIMESTEP: " << dt << endl;
  cout << "---------------------------" << endl;

  TimerList timers(1);
  timers.CreateTimer("Solve");

  timers.StartTimer("Solve");
  //compute analytic solution
  switch(mode)
    {
    case 0:
      if(dof > 1){
	cerr << "WARNING: Analytic solution only available for 1 DOF\n" << endl;
	return (-1);
      }
      cout << "Using Analytic Solution!" << endl;
      AnalyticSolution(dof, nsteps, dt, m, c, k, x, xd, xdd, &forces);
      break;
    case 1:
      cout << "Using RK 2nd Order Ralston!" << endl;
      RK2ndOrder(dof, nsteps, dt, m, c, k, x, xd, xdd, &forces, &bcs, 0);
      break;
    case 2:
      cout << "Using RK 2nd Order Trapezoidal!" << endl;
      RK2ndOrder(dof, nsteps, dt, m, c, k, x, xd, xdd, &forces, &bcs, 1);
      break;
    case 3:
      cout << "Using RK 3rd Order!" << endl;
      RK3rdOrder(dof, nsteps, dt, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 4:
      cout << "Using RK 4th Order!" << endl;
      RK4thOrder(dof, nsteps, dt, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 5:
      cout << "Using Classical Central Difference!" << endl;
      CentralDiff(dof, nsteps, dt, m, diagm, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 6:
      cout << "Using Half-step Central Difference!" << endl;
      HalfStepCentral(dof, nsteps, dt, m, diagm, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 7:
      cout << "Using Newmark-Beta - AVERAGE ACCELERATION!" << endl;
      gamma = 0.5;
      beta = 0.25;
      NewmarkBeta(dof, nsteps, dt, gamma, beta, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 8:
      cout << "Using Newmark-Beta - LINEAR ACCELERATION!" << endl;
      gamma = 0.5;
      beta = 1.0/6.0;
      NewmarkBeta(dof, nsteps, dt, gamma, beta, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 9:
      gamma = 0.5;
      beta = 1.0/12.0;
      cout << "Using Newmark-Beta - Fox-Goodwin!" << endl;
      NewmarkBeta(dof, nsteps, dt, gamma, beta, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    case 10:
      //theta > 1 with theta = 1.4 typical
      theta = 1.4;
      cout << "Using Wilson-theta!" << endl;
      WilsonTheta(dof, nsteps, dt, theta, m, c, k, x, xd, xdd, &forces, &bcs);
      break;
    default:
      cerr << "MODE: " << mode << " not valid -- exiting" << endl;
      break;
    }
  timers.StopTimer("Solve");
  timers.PrintAllTimers();

  //write solution
  WriteSolution(dof, nsteps, dt, x, xd, xdd);

  //write solution in vtk format
  double* temp = new double[param.nnodes*3];
  std::string* varNames = new std::string[3];
  varNames[0] = "PositionDx";
  varNames[1] = "PositionDy";
  varNames[2] = "PositionDz";
  std::string casename;

  //for visualizing movement
  double* rxyz = new double[param.nnodes*3];

#if 1
  for(t = 0; t < nsteps; t++){
    //load up solution array
    jj = 0;
    for(i = 0; i < param.nelem[0]; i++){
      int ndofpn = param.elements[i].GetNodeDOF();
      int ndof = param.elements[i].GetDxDOF();
      for(j = 0; j < param.elements[i].GetNnodes(); j++){
	int node = param.elements[i].nodes[j];
	int offset = param.nodeOffsetsDOF[node] + ndof;
	temp[jj*param.nnodes + node] = x[t*dof + offset];
      }
    }
    jj = 1;
    for(i = 0; i < param.nelem[0]; i++){
      int ndofpn = param.elements[i].GetNodeDOF();
      int ndof = param.elements[i].GetDyDOF();
      for(j = 0; j < param.elements[i].GetNnodes(); j++){
	int node = param.elements[i].nodes[j];
	int offset = param.nodeOffsetsDOF[node] + ndof;
	temp[jj*param.nnodes + node] = x[t*dof + offset];
      }
    }
    jj = 2;
    for(i = 0; i < param.nelem[0]; i++){
      int ndofpn = param.elements[i].GetNodeDOF();
      int ndof = param.elements[i].GetDzDOF();
      for(j = 0; j < param.elements[i].GetNnodes(); j++){
	int node = param.elements[i].nodes[j];
	int offset = param.nodeOffsetsDOF[node] + ndof;
	temp[jj*param.nnodes + node] = x[t*dof + offset];
      }
    }
    //load up displacement array for viz
    for(i = 0; i < param.nnodes; i++){
      rxyz[i*3 + 0] = temp[0*param.nnodes + i];
      rxyz[i*3 + 1] = temp[1*param.nnodes + i];
      rxyz[i*3 + 2] = temp[2*param.nnodes + i];
    }
    double scale = 1.0e2;
    param.MoveMesh(rxyz, rxyz, scale);

    ss.clear();
    ss.str("");
    ss << t;
    casename = "sol." + ss.str();
    WriteVTK_Ascii(casename, param.nelem, param.elements, param.nnodes, 
		   rxyz, temp, 3, varNames);
  }
#endif  

  delete [] m;
  delete [] diagm;
  delete [] c;
  delete [] k;
  delete [] x;
  delete [] xd;
  delete [] xdd;
  delete [] temp;

  return 0;
}


