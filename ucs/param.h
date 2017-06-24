#ifndef PARAM_H__
#define PARAM_H__

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>				
#include <sstream>
#include <vector>

#include "general.h"
#include "parameterParser.h"

//this is a list of the modes available to the solver
//these are mostly used for computational design
enum SolverMode
  {
    CFD,
    ObjectiveEval,
    Direct,
    Adjoint,
    CTSE,
    GridSmoothing,
    MeshSensitivity,
    FiniteDifference
  };

 
template <class Type>
class Param
{

public:
  Param();
  template <class Type2>
  Param(const Param<Type2>& toCopy);
  ~Param();

  //keep a copy of the casename -- used to look for restarts
  std::string spacename;
  std::string path;
  std::string casestring;
  std::string defaultChemDB;
  
  //mesh reordering to reduce bandwidth
  Bool reorder;

  //symmetry in directions x(0), y(1), z(2) -- manually zero the fluxes in these directions
  Int symmetry2D;

  //viscosity terms
  Int viscous;
  Int turbModel;
  Int turbModelSorder;

  //chemistry stuff
  Int chemModelId;
  Bool rxnOn;
  std::string chemDB;

  //reference values for non-dimensionalization
  //BASIC: 
  //-------------------------------
  Type ref_length;       //m
  Type ref_velocity;     //m/s
  Type ref_temperature;  //Kelvin
  Type ref_density;      //kg/m^3
  Type ref_viscosity;    // mu  (Pa. s) = (kg /(m. s)
  Type ref_k;            // W/(m.K)
  //DERIVED:
  //--------------------------------
  Type ref_time;         //s
  Type ref_enthalpy;     //(m^2/s^2)/Ec
  Type ref_pressure;     //(kg /( m. s^2))

  Type Re;            // Reynolds number - rho*U*L/Mu
  Type Pr;            // Prandtl number - Cp*Mu/k
  Type PrT;           // Turbulent prandtl number ~0.7-0.9

  Bool scaleMesh;     // If set true, we scale the mesh by ref_length before computation

  Type gamma;         // ratio of specific heats
  Bool dynamicCFL;    // turns on CFL dynamic scaling
  Type cflStart;      // CFL number to start ramping from
  Type cflPrev;       // CFL from previous timestep
  Type cfl;           // CFL number
  Int cflRampingSteps;// number of steps to ramp used cfl to "cfl" variable
  Bool useLocalTimeStepping; // 0 - false, 1 - true
  Bool enableVNN;     // von neumann diffusive stability condition - useful for very low Reynolds numbers
  Type VNN;           // von neumann diffusive stability multiplier

  Type time;            // total elapsed time
  Type dt;              // minimum time step (unsteady only); if negative, use
                        // local time stepping
  Bool pseudotimestepping; // use pseudotimestepping
  Type dtau;            // minimum pseudo time step
  int torder;           // temporal order (1 or 2)
  int sorder;           // spatial order
  Type chi;             // 0.5 - quadratic reconstruction, 0.0 - linear reconstruction
  int nFirstOrderSteps; // number of first order spatial steps to take before switching
                        // to a higher order scheme
  int limiter;          // 0 - no limiter
                        // 1 - Barth
                        // 2 - Venkatakrishnan
                        // 3 - Modified Venkatakrishnan
  int limiterRefresh;   // number of steps to wait to refresh the limiter
                        // 1 - immediately refresh
                        // 2 - every other step 
                        // etc.
  int limiterFreeze;    // step to freez limiter updates at

  int jacobianFreq;         // frequency to update Jacobians
  int boundaryJacEval;      // 0 - exact, 1 - one-sided (approx.)
  int boundaryJacType;      // 0 - one sided, 1 - central diff, 2 - complex
  int fieldJacType;         // 0 - one sided, 1 - central diff, 2 - complex
  Bool hojac;                // 0 - first order, 1 - order of spatial solver

  //flag which tells us if we require a complex
  //eqnset and chem model to be built for jacobians
  Bool requireComplex;

  int nSgs;                    // number of sgs sweeps if 0 use explicit method


  Bool useRestart;       // if true look for restart file (default 0)
  int writeRestartStep;  // number of steps to wait before writing new 
                         // restart file (i.e. 100 writes every 100 steps)
  Bool preserveRestartCounters; //if true, preserve the iteration counters
                                //just like a restart never happened

  Bool solutionTagStep;    //controls number of solution files to keep around
  int solutionWrite;       //number of steps to wait before writing solution files 
  
  int customIcId;       // >0 use ic defined in customics.cpp

  std::vector<Type> flowdir;      //flow direction vector
  std::vector<Type> liftdir;      //lift direction vector
  std::vector<Type> dragdir;      //drag direction vector

  std::vector<Type> gravdir;      //gravity direction vector
  Type gravity;                   //gravity magnitude
  Bool gravity_on;                //gravity switch

  std::vector<Type> massfractions;   //mass fractions of species for reacting flows
  std::vector<Type> molefractions;   //mole fractions of species for reacting flows (optional setting)

  Int eqnset_id;     //see eqnset_defines.h for enumeration of types

  Bool no_cvbc;

  Type beta;         // artificial compressibility parameter (for incompressible)
  Type betaMin;      // minimum preconditioner value (otherwise set to Ma^2 in most cases)
  Int flux_id;       //flux type id -- see eqnset.h
                     //0 - Roe
                     //1 - HLLC

  Int gradType;  //0 - LSQ
                 //1 - Green Gauss

  //movement routine to use
  //0 - none
  //1 - pitching solid body rotation
  //2 - pitching dynamic mesh
  Int movement;   
  //enable the Geometric Conservation Law for changing grids (GCL)
  Int gcl;

  //This section is Computational Design stuff only
  Int mode;                   //keeps track of the mode the solver is running in
  Int designSolver;           //0 - SGS
                              //1 - GMRES
  Int designPrecond;          //0 - none
  Int designNsgs;             //number of sgs iterations
  Int designRestarts;         //GMRES - number of restarts 
  Int designSearchDir;        //search directions before GMRES restart
  Int objFuncId;              //id of the objective function
  std::vector<Type> sensTarget;  //list of target values for the sensors defined in order
  Int sensEqn;                //number of equation (zero based) which we give sensor targets for

  //list of field names which are written to solution file
  std::vector<std::string> fieldsRequested;

  //section related to gaussian source terms
  Bool gaussianSource;
  Int gaussianEqn;
  Int gaussianBCid;
  Type gaussianAmpl;
  Type gaussianVelAmpl;
  Type gaussianXloc;
  Type gaussianYloc;

  //section related to error transport equations
  Bool errorTransport;

  Type velocityStart;       //velocity to start ramping from
  Type velocity;            //velocity in farfield
  Int velocityRampingSteps; //number of steps used to ramp velocity to "velocity" variable

  void SetupParams();

  Int ReadSpace(std::string spaceName, std::ifstream& fin, std::streampos locBegin, std::streampos locEnd);
  //compute implicitly declared parameters
  void PostCompute();
  void PrintSolverParams();
  void PrintAllParams();

  void UpdateCFL(Int iter, Type deltaResidual);
  Type GetCFL();
  Type GetVelocity(Int iter);

  void UpdateForDesign();
  
  const ParameterString& GetParameterString(std::string name) const;
  const ParameterStringList& GetParameterStringList(std::string name) const;
  const ParameterBool& GetParameterBool(std::string name) const;
  const Parameter<Int>& GetParameterInt(std::string name) const;
  const Parameter<Type>& GetParameterType(std::string name) const;
  const ParameterList<Type>& GetParameterListType(std::string name) const;
  const ParameterEnum& GetParameterEnum(std::string name) const;

 protected:

 private:

  Int ParseLine(std::string& line);

  //a list of all of the possible parameters which are checked for
  std::vector<ParameterBool> paramListBool;
  std::vector<Parameter<Int> > paramListInt;
  std::vector<Parameter<Type> > paramListReal;
  std::vector<ParameterList<Type> > paramListRealList;
  std::vector<ParameterEnum> paramListEnum;
  std::vector<ParameterString> paramListString;
  std::vector<ParameterStringList> paramListStringList;
};


//This is the master param file reading call, will create and load a vector
//of param files which correspond to each solution space defined in the param file
template <class Type> 
Int ReadParamFile(std::vector<Param<Type>*>& paramFileList, std::string casename, std::string pathname, 
		  Bool verbose = true);

//include implementations
#include "param.tcc"

#endif
