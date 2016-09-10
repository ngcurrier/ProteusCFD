#include "fluid_structure.h"
#include "geometry.h"
#include "h5layer.h"

namespace STRUCTDYN{

TypicalSection::TypicalSection(Param<double>* param, std::string name, TemporalControl<double>& temporalControl) :
  SolutionSpaceBase<double>(name, temporalControl)
{
  this->param = param;
  //use average acceleration O(dt^2)
  gamma = 0.5;
  beta = 0.25;
  dof = 2;

  m = new double[dof*dof];
  k = new double[dof*dof];
  c = new double[dof*dof];
  rhs = new double[dof];

  dx = new double[dof];
  x_n = new double[dof];
  xd_n = new double[dof];
  xdd_n = new double[dof];

  this->AddScalarField("Plunge");
  this->AddScalarField("Pitch");

  return;
}

TypicalSection::~TypicalSection()
{
  delete [] m;
  delete [] k;
  delete [] c;
  delete [] rhs;
  delete [] dx;
  delete [] x_n;
  delete [] xd_n;
  delete [] xdd_n;

  return;
}

void TypicalSection::PreIterate()
{
  int i;

  //3 - PAPA flutter case - NACA0012
  //6 - Isogai flutter case - NACA64A010
  if(param->movement == 3){
    //this is for the PAPA mechanism flutter case - NACA0012
    ma = 87.067225;  //conversion slugs -> kg multiply by 14.59
    //mass should be the mass per unit span, the PAPA case is 32 in. span = 0.8128 m
    ma = ma / 0.8128;
    chord = param->ref_length;
    kalpha = 3927.96; //conversion lbs*ft/rad -> N*m/rad multiply by 1.356
    kh = 38806.732;  //conversion lbs/ft -> N/m multiply by 14.59
    xalpha = 0.0;
    //icg = iea - m(b*x_alpha)^2
    //iea = k_alpha/omega_alpha^2 -> omega_alpha is natural frequency (rot.) at zero airspeed
    icg = 3.679;
    rhoinf = param->ref_density;
    double V = param->GetVelocity(this->iter)*param->ref_velocity;
    double Mach = V;
    uinf = Mach;
  }
  else if(param->movement == 6){
    //this is for the Isogai flutter case - NACA64A010
    ma = 57.726765;
    chord = param->ref_length; //should be set to 1.0 meters
    kalpha = 502222.855;
    kh = 577267.65;
    xalpha = 1.8;  //(static mass moment per unit span )/ (mass * b)
    icg = 3.46360595;
    //rhoinf should be 1.225 kg/m^3
    rhoinf = param->ref_density; 
    double V = param->GetVelocity(this->iter)*param->ref_velocity;
    double Mach = V;
    uinf = Mach;
  }
  else{
    std::cerr << "Invalid Typical Section case type" << std::endl;
    std::cerr << "Set movement == 3 for PAPA and movement == 6 for Isogai setups" << std::endl;
    std::cerr << "movement flag found is: " << param->movement << std::endl;
    return;
  }

  //semi-chord
  this->b = chord/2.0;
  //mass moment of inertia about elastic axis (parallel axis thm.)
  this->iea = icg + ma*(b*xalpha)*(b*xalpha);
  //natural frequency rotation
  this->omega_alpha = sqrt(kalpha/iea);
  //nondimensionalize the time by natural frequency of rotation
  //the dt parameter is assumed non-dimensional time, so here we multiply it by the reference 
  //quantity to get dimensional time before non-dimensionalizing again
  this->dt = (param->dt*param->ref_time)*this->omega_alpha;

  BuildSystem();

  if(!param->useRestart){
    this->iter = 1;
    //if not using restarting, use standard ICs
    int dof = 2;
    double* x = new double[dof];
    double* xd = new double[dof];
    for(int kk = 0; kk < dof; kk++){
      x[kk] = xd[kk] = 0.0;
    }
    if(param->movement == 3){
      //perturb the airfoil for the PAPA case, we force Isogai through a few cycles
      //this is done elsewhere
      xd[0] = xd[1] = 0.0002;
    }
    IC(x, xd);
    delete [] x;
    delete [] xd;
  }
  else{
    ReadRestartFile();
  }


  std::cout << "Computed uncoupled natural frequencies: Omega_alpha - " 
	    << this->omega_alpha << " Omega_h - " << sqrt(kh/ma) << std::endl;
  std::cout << "Nondimensional structural dt: " << this->dt << std::endl;
  std::cout << "Dimensional structural dt: " << param->dt*param->ref_time << std::endl;
  std::cout << "Moment of inertia - elastic axis: " << this->iea << std::endl;
  std::cout << "Radius of gyration: " << sqrt(iea/(ma*b*b)) << std::endl;
  double pi = acos(-1.0);
  double mu = ma/(pi*rhoinf*b*b);
  std::cout << "Airfoil mass ratio: " << mu << std::endl;
  std::cout << "Flutter speed index: " << uinf/(b*omega_alpha*sqrt(mu)) << std::endl;

  return;
}

void TypicalSection::IC(double* x, double* xd)
{
  Int i, j;
  for(i = 0; i < dof; i++){
    dx[i] = x_n[i] = xd_n[i] = 0.0;
  }

  for(i = 0; i < dof; i++){
    x_n[i] = x[i];
    xd_n[i] = xd[i];
  }

  double* mtemp = new double[dof*dof];
  double* t1 = new double[dof];
  double* t2 = new double[dof];
  double* diag = new double[dof];

  //build a temp rhs with cl and cm equal to zero
  ScalarField<double>& clf = GetScalarField("CL");
  ScalarField<double>& cmf = GetScalarField("CM");
  clf.SetField(0.0);
  cmf.SetField(0.0);
  BuildRHS();

  //initialize acceleration for first time step
  MatVecMult(c, xd_n, t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] = rhs[j] - t2[j];
  }

  MatVecMult(k, x_n, t2, dof);
  for(j = 0; j < dof; j++){
    t1[j] -= t2[j];
  }

  memcpy(mtemp, m, sizeof(double)*dof*dof);

  //compute cholesky factorization
  Cholesky(mtemp, diag, dof);
  CholeskySolve(mtemp, diag, xdd_n, t1, dof);

  //this is kinda HACKY... watch out...
  //b/c of the way acclerations are computed above.. this method blows
  //up for higher than avg. accelerations (beta/gamma)
  //i.e. the velocity and acceleration terms for static ic's go nuclear
  //zero them if initial velocity is zero
  for(j = 0; j < dof; j++){
    if(xd_n[j] == 0.0 && this->iter == 1){
      xdd_n[j] = 0.0;
    }
  }

  std::cout << "INITIAL CONDITIONS: STRUCTURAL COUPLING" << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  for(i = 0; i < dof; i++){
    std::cout << x_n[i] << " " << xd_n[i] << " " << xdd_n[i] << std::endl;
  }
  std::cout << std::endl;

  delete [] mtemp;
  delete [] t1;
  delete [] t2;
  delete [] diag;

  return;
}

void TypicalSection::WriteRestartFile()
{
  std::string filename = param->path+param->spacename + ".0." + "rs";
  
  std::cout << "RESTART WRITEFILE: writing restart file --> " << filename << std::endl;

  hid_t h5out = HDF_OpenFile(filename, 1);
  HDF_WriteScalar(h5out, "/Solver State/", "Iteration Count", &this->iter);
  HDF_WriteArray(h5out, "/Fields/", "dx", dx, dof);
  HDF_WriteArray(h5out, "/Fields/", "x_n", x_n, dof);
  HDF_WriteArray(h5out, "/Fields/", "xd_n", xd_n, dof);
  HDF_WriteArray(h5out, "/Fields/", "xdd_n", xdd_n, dof);
  HDF_CloseFile(h5out);

  return;
}

void TypicalSection::ReadRestartFile()
{
  std::string filename =  param->path+param->spacename + ".0." + "rs";
  std::cout << "RESTART READFILE: reading restart file --> " << filename << std::endl;
  
  hid_t h5in = HDF_OpenFile(filename, 0);
  HDF_ReadScalar(h5in, "/Solver State/", "Iteration Count", &this->iter);
  HDF_ReadArray(h5in, "/Fields/", "dx", &dx, &dof);
  HDF_ReadArray(h5in, "/Fields/", "x_n", &x_n, &dof);
  HDF_ReadArray(h5in, "/Fields/", "xd_n", &xd_n, &dof);
  HDF_ReadArray(h5in, "/Fields/", "xdd_n", &xdd_n, &dof);
  HDF_CloseFile(h5in);

  //echo the input
  std::cout << "TYPICAL SECTION: Restarting dx: " << dx[0] << " " << dx[1] << std::endl;
  std::cout << "TYPICAL SECTION: Restarting pos: " << x_n[0] << " " << x_n[1] << std::endl;
  std::cout << "TYPICAL SECTION: Restarting velocity: " << xd_n[0] << " " << xd_n[1] << std::endl;
  std::cout << "TYPICAL SECTION: Restarting acceleration: " << xdd_n[0] << " " << xdd_n[1] << std::endl;

  return;
}

void TypicalSection::BuildSystem() 
{
  //natural frequency plunge
  double omega_h = sqrt(kh/ma);
  //radius of gyration
  double ralpha = sqrt(iea/(ma*b*b));
  double temp;

  int n = dof;

  for(int i = 0; i < dof*dof; i++){
    c[i] = 0.0;
  }

  //first row m
  m[0*n + 0] = 1.0;
  m[0*n + 1] = xalpha;
  //second row m
  m[1*n + 0] = xalpha;
  m[1*n + 1] = ralpha*ralpha;

  //first row k
  temp = omega_h/omega_alpha;
  k[0*n + 0] = temp*temp; 
  k[0*n + 1] = 0.0;
  //second row k
  k[1*n + 0] = 0.0;
  k[1*n + 1] = ralpha*ralpha;

  return;
}

void TypicalSection::BuildRHS()
{
  ScalarField<double>& clf = GetScalarField("CL");
  ScalarField<double>& cmf = GetScalarField("CM");

  double cl = clf.GetField();
  double cm = cmf.GetField();

  //build rhs
  double pi = acos(-1.0);
  //airfoil mass ratio
  double mu = ma/(pi*rhoinf*b*b);
  //reduced frequency
  double kc = omega_alpha*chord/uinf;
  double temp = 4.0/(pi*mu*kc*kc);
 
  rhs[0] = temp*(-cl);
  rhs[1] = temp*(2.0*cm);

  std::cout << "TYPICAL_SECTION: CL - " << cl << std::endl;
  std::cout << "TYPICAL_SECTION: CM - " << cm << std::endl;

  return;
}

void TypicalSection::NewtonIterate()
{
  ScalarField<double>& pitch = GetScalarField("Pitch");
  ScalarField<double>& plunge = GetScalarField("Plunge");

  if(iter < 50){
    plunge.SetField(0.0);
    pitch.SetField(0.0);
    return;
  }

  BuildRHS();

  //this gives us two driven cycles
  double ndt = (double)iter*dt;
  double x_np1[dof];
  if(ndt <= 12.5){
    ForcedMovement(dx, x_n, x_np1, xd_n, xdd_n);
  }
  else{
    //the unknowns in the non-dimensional system are h/b-dof#1 and alpha-dof#2
    NewmarkBetaDx(dof, dt, gamma, beta, m, c, k, dx, x_n, xd_n, xdd_n, rhs);
  }

  double h = (x_n[0] + dx[0])*b;
  double alpha = (x_n[1] + dx[1]);

  plunge.SetField(h);
  pitch.SetField(alpha);

  std::cout << "TYPICAL_SECTION: alpha - " << alpha << std::endl;
  std::cout << "TYPICAL_SECITON: plunge - " << h << std::endl;
		
  return;
}

void TypicalSection::PostTimeAdvance()
{
  int i;
  double x_np1[dof];
  double xd_np1[dof];
  double xdd_np1[dof];

  NewmarkBetaUpdateVelAcc(dof, dt, gamma, beta, dx, xd_n, xdd_n, xd_np1, xdd_np1);

  //this gives us four driven cycles
  double ndt = (double)iter*dt;
  if(ndt <= 12.5){
    ForcedMovement(dx, x_n, x_np1, xd_np1, xdd_np1);
  }

  //copy down the new solution to the old solution location
  for(i = 0; i < dof; i++){
    x_n[i] = x_n[i] + dx[i];
    xd_n[i] = xd_np1[i];
    xdd_n[i] = xdd_np1[i];
  }
  iter++;

  
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

  return;
}

void TypicalSection::Print(std::ostream& sout)
{
  sout << iter << ": " << iter*dt << " " << x_n[0] << " " << x_n[1] << " " 
       << xd_n[0] << " " << xd_n[1] << " " << xdd_n[0] 
       << " " << xdd_n[1] << std::endl;
  return;
}


void TypicalSection::ForcedMovement(double* dx, double* x_old, double* x, double* xd, double* xdd)
{
  double ndt = (double)iter*dt;
  double alpha = 0.017453292*sin(ndt);
  double dalpha = 0.017453292*cos(ndt);
  double ddalpha = -0.017453292*sin(ndt);

  //plunge is frozen
  dx[0] = 0.0;
  x[0] = 0.0;
  xd[0] = 0.0;
  xdd[0] = 0.0;

  //pitch is forced
  x[1] = alpha;
  dx[1] = alpha - x_old[1];
  xd[1] = dalpha;
  xdd[1] = ddalpha;

  return;
}

Structure::Structure(Param<double>* param, std::string name, TemporalControl<double>& temporalControl) :
  SolutionSpaceBase<double>(name, temporalControl)
{
  this->param = param;
  connect_alloc = false;
  rays = NULL;
  feaxyz = NULL;
  elem = NULL;

  return;
}

Structure::~Structure()
{
  if(connect_alloc){
    delete [] rays;
    delete [] elem;
    delete [] feaxyz;
  }
  return;
}


void Structure::Init()
{
  //TODO: fill this in 
  return;
}


void Structure::BuildFieldConnectivity(double* xyz, int npts)
{
  int i, j;

  nwettedpts = npts;
  //do some memory allocation
  if(connect_alloc){
    std::cerr << "WARNING: Memory already allocated BuildFieldConnectivity()!! ERROR!" 
	      << std::endl;
    return;
  }
  rays = new double[npts*3];
  elem = new int[npts];
  feaxyz = new double[npts*3];
  connect_alloc = true;


  for(i = 0; i < npts; i++){
    elem[i] = GetClosestElement(&xyz[i*3], &this->rays[i*3]);
  }
  
  //get the point we are interpolating to and store it
  for(i = 0; i < npts; i++){
    for(j = 0; j < 3; j++){
      this->feaxyz[i*3 + j] = xyz[i*3 + j] - this->rays[i*3 + j];
    }
  }
  
  return;
}

void Structure::GetPointsDisplacements(double* dxyz)
{
  //The idea here is to create a ray from the surface (CFD) point to the nearest point
  //on any element and then to use it's rotation and translation (FE sol'n) to move
  //the surface node in a conservative way

  int i, j;

  int maxDOF = 6;
  int ndof = 0;
  double values[maxDOF];

  for(i = 0; i < nwettedpts; i++){
    
    sparam->elements[this->elem[i]].InterpolatePoint(&this->feaxyz[i*3], sparam, x_np1, 
						    &ndof, values);
    if(ndof != 6){
      std::cout << "WARNING: Feeling kinda funny about not having a 6 DOF element--- HELP ME!"
		<< std::endl;
    }
    //compute displacements, this assumes the DOF stored in values are in the 
    //following order Dx, Dy, Dz, Rx, Ry, Rz
    double r = Magnitude(&this->rays[i*3]);
    //the storage pattern below is an assumption, might be wrong in some non-6-DOF elements
    double dx = values[0];
    double dy = values[1];
    double dz = values[2];
    double rx = values[3];
    double ry = values[4];
    double rz = values[5];

    dxyz[i*3 + 0] = dx + r*cos(rz) - r + r*sin(ry);
    dxyz[i*3 + 1] = dy + r*cos(rx) - r + r*sin(rz);
    dxyz[i*3 + 2] = dz + r*cos(ry) - r + r*sin(rx);
  }

  return;
}

int Structure::GetClosestElement(double* xyz, double* radVector)
{
  int i, j, k;
  int elem;

  //coordinate which we assume axis is along
  //this will have to change if we use anymore than 1D beams in the future
  int fixed_dir = 1;  //y - direction
  
  for(j = 0; j < sparam->nelem[0]; j++){
    int nnodes = sparam->elements[j].GetNnodes();
    int n1 = sparam->elements[j].nodes[0];
    int n2 = sparam->elements[j].nodes[1];
    //check to see if node is within the bounds of the element
    if((sparam->xyz[n1*3 + fixed_dir] < xyz[fixed_dir] && 
	sparam->xyz[n2*3 + fixed_dir] > xyz[fixed_dir]) ||
       (sparam->xyz[n1*3 + fixed_dir] > xyz[fixed_dir] && 
	sparam->xyz[n2*3 + fixed_dir] < xyz[fixed_dir])){
      elem = j;
      double dist = DistanceLinePoint(&sparam->xyz[n1*3], &sparam->xyz[n2*3], xyz, radVector);
      break;
    }
  }

  return elem;
}

void Structure::ApplyCFDForce(double* force)
{
  int i;
  int elemid;
  int maxDOF = 6;
  double feaForces[maxDOF];

  for(i = 0; i < nwettedpts; i++){
    //get containing element
    elemid = elem[i];
    //check to make sure we have a full 6DOF element
    int ndof = this->sparam->elements[elemid].GetNodeDOF();
    if(ndof != 6){
      std::cout << "WARNING: Feeling kinda funny about not having a 6 DOF element--- HELP ME!"
		<< std::endl;
    }
    //compute forces, this assumes the DOF stored in values are in the 
    //following order Fx, Fy, Fz, Mx, My, Mz
    feaForces[0] = force[i*3 + 0];
    feaForces[1] = force[i*3 + 1];
    feaForces[2] = force[i*3 + 2];
    //compute moments
    //M = r X F
    double Moment[3];
    CrossProduct(&rays[i*3], &force[i*3], &feaForces[3]);
    
    //interpolate the forces/moments from the element to 
    //the element nodes where they can actually be applied
    //this call also accumulates the forces to the correct 
    //DOF in the rhs
    this->sparam->elements[elemid].InterpolateValuesToNodes(&feaxyz[i*3], sparam, feaForces, 
							   ndof, this->rhs);
  }

  return;
}

void Structure::ComputeDx()
{
  //the unknowns in the non-dimensional system are h/b-dof#1 and alpha-dof#2
  NewmarkBetaDx(dof, dt, gamma, beta, m, c, k, dx, x_n, xd_n, xdd_n, rhs);
		
  return;
}

void Structure::WriteRestart(std::ofstream& fout)
{
  return;
}

void Structure::ReadRestart(std::ifstream& fin)
{
  return;
}

void Structure::NextTimestep()
{
  int i;
  double x_np1[dof];
  double xd_np1[dof];
  double xdd_np1[dof];

  NewmarkBetaUpdateVelAcc(dof, dt, gamma, beta, dx, xd_n, xdd_n, xd_np1, xdd_np1);

  //copy down the new solution to the old solution location
  for(i = 0; i < dof; i++){
    x_n[i] = x_n[i] + dx[i];
    xd_n[i] = xd_np1[i];
    xdd_n[i] = xdd_np1[i];
  }
  this->iter++;

  return;
}

}
