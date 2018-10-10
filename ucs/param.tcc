
#include "eqnset_defines.h"
#include "strings_util.h"
#include "vars_id.h"
#include "exceptions.h"
#include "portFileio.h"
#include "macros.h"
#include <sstream>
#include <stdlib.h>

#define REF_PRESSURE_DEFAULT 1.0
#define REF_DENSITY_DEFAULT 1.0
#define R_UNIV 8.314459848      // J/mol.K


template <class Type>
Param<Type>::Param()
{
  SetupParams();
}

template <class Type> template <class Type2>
Param<Type>::Param(const Param<Type2>& toCopy)
{
  //call to insure internals are ready
  SetupParams();
  path = toCopy.path;
  spacename = toCopy.spacename;
  casestring = toCopy.casestring;
  mode = toCopy.mode;

  for(std::vector<ParameterString>::iterator it = paramListString.begin();
      it != paramListString.end(); ++it){
    ParameterString& p = *it;
    const ParameterString& ptc = toCopy.GetParameterString(p.GetKey());
    p.Copy(ptc);
  }
  for(std::vector<ParameterStringList>::iterator it = paramListStringList.begin();
      it != paramListStringList.end(); ++it){
    ParameterStringList& p = *it;
    const ParameterStringList& ptc = toCopy.GetParameterStringList(p.GetKey());
    p.Copy(ptc);
  }
  for(std::vector<ParameterBool>::iterator it = paramListBool.begin();
      it != paramListBool.end(); ++it){
    ParameterBool& p = *it;
    const ParameterBool& ptc = toCopy.GetParameterBool(p.GetKey());
    p.Copy(ptc);
  }
  for(std::vector<Parameter<Int> >::iterator it = paramListInt.begin();
      it != paramListInt.end(); ++it){
    Parameter<Int>& p = *it;
    const Parameter<Int>& ptc = toCopy.GetParameterInt(p.GetKey());
    p.Copy(ptc);
  }
  for(typename std::vector<Parameter<Type> >::iterator it = paramListReal.begin();
      it != paramListReal.end(); ++it){
    Parameter<Type>& p = *it;
    const Parameter<Type2>& ptc = toCopy.GetParameterType(p.GetKey());
    p.Copy(ptc);
  }
  for(typename std::vector<ParameterList<Type> >::iterator it = paramListRealList.begin();
      it != paramListRealList.end(); ++it){
    ParameterList<Type>& p = *it;
    const ParameterList<Type2>& ptc = toCopy.GetParameterListType(p.GetKey());
    p.Copy(ptc);
  }
  for(std::vector<ParameterEnum >::iterator it = paramListEnum.begin();
      it != paramListEnum.end(); ++it){
    ParameterEnum& p = *it;
    const ParameterEnum& ptc = toCopy.GetParameterEnum(p.GetKey());
    p.Copy(ptc);
  }

  //we call this b/c some post-processing is done in here
  PostCompute();
}

template <class Type>
Param<Type>::~Param()
{ }

template <class Type>
void Param<Type>::SetupParams()
{
  //These are internal non-user parameters, used for runtime monitoring, etc.
  this->chi = 0.0;
  this->requireComplex = false;
  this->PrT = 0.85;  //average value??? no idea if this is okay varies from 0.7-0.9
  this->ref_time = 1.0;
  this->ref_enthalpy = 1.0;
  this->ref_density = 1.0;
  this->ref_pressure = 1.0;
  this->chemModelId = 0;   //model 0 is defined from file
  this->Re = 1.0;
  this->viscous = 0;

  const char* homeDir = getenv("HOME");
  std::string homeStr(homeDir);
  defaultChemDB = homeStr + "/.proteusCFD/database/chemdb.hdf5";

  //SET THE STRING VALUED PARAMETERS
  paramListString.push_back(ParameterString("chemicalDatabase", &chemDB, defaultChemDB));

  //SET THE REAL VALUED LIST PARAMETERS
  std::vector<Type> val(3);
  //inflow in x-direction default
  val[0] = 1.0;
  val[1] = 0.0;
  val[2] = 0.0;
  paramListRealList.push_back(ParameterList<Type>("flowDirection", &flowdir, val, -1.0, 1.0, 3, 3));
  //drag in x-direction default
  paramListRealList.push_back(ParameterList<Type>("dragDirection", &dragdir, val, -1.0, 1.0, 3, 3));
  //lift in y-direction default
  val[0] = 0.0;
  val[1] = 1.0;
  val[2] = 0.0;
  paramListRealList.push_back(ParameterList<Type>("liftDirection", &liftdir, val, -1.0, 1.0, 3, 3));
  //gravity in negative z-direction default
  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = -1.0;
  paramListRealList.push_back(ParameterList<Type>("gravityDirection", &gravdir, val, -1.0, 1.0, 3, 3));
  val.clear();
  val.resize(0);
  paramListRealList.push_back(ParameterList<Type>("massFractions", &massfractions, val, 0.0, 1.0, 999, 1));
  paramListRealList.push_back(ParameterList<Type>("moleFractions", &molefractions, val, 0.0, 1.0, 999, 1));
  paramListRealList.push_back(ParameterList<Type>("sensorTarget", &sensTarget, val, 0.0, 9999999999.0, 999, 1));

  //SET THE STRING LIST VALUED PARAMETERS
  paramListStringList.push_back(ParameterStringList("fieldsRequested", &fieldsRequested, "variableQ"));

  //SET THE ENUM VALUED PARAMETERS
  std::vector<std::string> eqnList;
  for(Int i = 0; i < NUM_EQN_SETS; i++){
    //eqnSet and fluxType strings are included from eqnset.h
    //look there to change defintion of those keywords
    eqnList.push_back(eqnSets[i]);
  }
  paramListEnum.push_back(ParameterEnum("equationSet", eqnList, &this->eqnset_id, CompressibleEuler));
  std::vector<std::string> fluxList;
  for(Int i = 0; i < NUM_FLUX_TYPES; i++){
    fluxList.push_back(fluxTypes[i]);
  }
  paramListEnum.push_back(ParameterEnum("fluxType", fluxList, &this->flux_id, RoeFluxType));

  //SET THE BOOLEAN VALUED PARAMTERS
  paramListBool.push_back(ParameterBool("localTimeStepping", &this->useLocalTimeStepping, true));
  paramListBool.push_back(ParameterBool("higherOrderJacobians", &this->hojac, false));
  paramListBool.push_back(ParameterBool("noCVBC", &this->no_cvbc, false));   //default to using CVBCs where available
  paramListBool.push_back(ParameterBool("reorderMesh", &this->reorder, true));  //default to doing reordering
  paramListBool.push_back(ParameterBool("useRestart", &this->useRestart, false));
  paramListBool.push_back(ParameterBool("preserveRestartCounters", &this->preserveRestartCounters, true));
  paramListBool.push_back(ParameterBool("reactionsOn", &this->rxnOn, true));
  paramListBool.push_back(ParameterBool("solutionTagStep", &this->solutionTagStep, false));
  paramListBool.push_back(ParameterBool("gaussianSourceOn", &this->gaussianSource, false));
  paramListBool.push_back(ParameterBool("scaleMesh", &this->scaleMesh, false));
  paramListBool.push_back(ParameterBool("errorTransport", &this->errorTransport, false));
  paramListBool.push_back(ParameterBool("dynamicCFL", &this->dynamicCFL, false));
  paramListBool.push_back(ParameterBool("enableVNN", &this->enableVNN, false));
  paramListBool.push_back(ParameterBool("enableGravity", &gravity_on, false));

  //SET THE INTEGER VALUED PARAMETERS
  paramListInt.push_back(Parameter<Int>("turbulenceModel", &this->turbModel, 0, 0, 999)); //default to laminar
  paramListInt.push_back(Parameter<Int>("turbulenceModelSpatialOrder", &this->turbModelSorder, 0, 0, 3)); 
  paramListInt.push_back(Parameter<Int>("numberSGS", &this->nSgs, 5, 0, 999));
  paramListInt.push_back(Parameter<Int>("timeOrder", &this->torder, 0, 0, 2)); 
  paramListInt.push_back(Parameter<Int>("spatialOrder", &this->sorder, 0, 0, 3));
  paramListInt.push_back(Parameter<Int>("firstOrderSteps", &this->nFirstOrderSteps, 0, 0, 999999));
  paramListInt.push_back(Parameter<Int>("limiterRefresh", &this->limiterRefresh, 1, 1, 999999));
  paramListInt.push_back(Parameter<Int>("limiterFreeze", &this->limiterFreeze, 999999, 0, 999999));
  paramListInt.push_back(Parameter<Int>("limiter", &this->limiter, 0, 0, 5));
  paramListInt.push_back(Parameter<Int>("jacobianUpdateFrequency", &this->jacobianFreq, 1, 1, 999999));
  paramListInt.push_back(Parameter<Int>("jacobianBoundaryEval", &this->boundaryJacEval, 0, 0, 2));
  paramListInt.push_back(Parameter<Int>("jacobianBoundaryType", &this->boundaryJacType, 0, 0, 2));
  paramListInt.push_back(Parameter<Int>("jacobianFieldType", &this->fieldJacType, 0, 0, 2));
  paramListInt.push_back(Parameter<Int>("gcl", &this->gcl, 0, 0, 2));
  paramListInt.push_back(Parameter<Int>("rampVelocity", &this->velocityRampingSteps, 0, 0, 999999));
  paramListInt.push_back(Parameter<Int>("movement", &this->movement, 0, 0, 999));  //default to no movement
  paramListInt.push_back(Parameter<Int>("restartWrite", &this->writeRestartStep, 0, 0, 999999));
  paramListInt.push_back(Parameter<Int>("customIcId", &this->customIcId, 0, 0, 99));
  paramListInt.push_back(Parameter<Int>("solutionWrite" , &this->solutionWrite, 0, 0, 9999999)); 
  paramListInt.push_back(Parameter<Int>("gradientType", &this->gradType, 0, 0, 1));
  paramListInt.push_back(Parameter<Int>("rampCFL", &this->cflRampingSteps, 0, 0, 999999));
  paramListInt.push_back(Parameter<Int>("symmetry2D", &this->symmetry2D, -1, 0, 3));
  paramListInt.push_back(Parameter<Int>("gaussianEqn", &this->gaussianEqn, 0, 0, 999));
  paramListInt.push_back(Parameter<Int>("gaussianBCid", &this->gaussianBCid, 0, 0, 999));

  //SET THE REAL VALUED PARAMETERS
  paramListReal.push_back(Parameter<Type>("startingCFL", &this->cflStart, 0.5, 0.0, 10000.0));
  paramListReal.push_back(Parameter<Type>("CFL", &this->cfl, 5.0, 0.0, 10000.0));
  paramListReal.push_back(Parameter<Type>("VNN", &this->VNN, 20.0, 0.0, 10000.0));
  paramListReal.push_back(Parameter<Type>("timeStep", &this->dt, -1.0, -1.0, 10000.0));
  paramListReal.push_back(Parameter<Type>("velocity", &this->velocity, 1.0, 0.0, 999.0));
  paramListReal.push_back(Parameter<Type>("startingVelocity", &this->velocityStart, 0.0, 0.0, 999.0));
  paramListReal.push_back(Parameter<Type>("MW", &this->MW, 28.966, 0.0, 999.0));  // default to air
  paramListReal.push_back(Parameter<Type>("Cp", &this->Cp, 1006.43, 0.0, 3000.0)); // default to air
  paramListReal.push_back(Parameter<Type>("betaMin", &this->betaMin, 0.0, 0.0, 999.0));
  paramListReal.push_back(Parameter<Type>("beta", &this->beta, 15.0, 0.0, 999.0)); 
  paramListReal.push_back(Parameter<Type>("prandtlNumber", &this->Pr, 0.72, 0.0, 999.0));
  paramListReal.push_back(Parameter<Type>("thermalConductivity", &this->kThermalConductivity, 1.0, 0, 99999.0));
  paramListReal.push_back(Parameter<Type>("specificHeat", &this->cpSpecificHeat, 1.0, 0, 9999.0)); 
  paramListReal.push_back(Parameter<Type>("density", &this->rhoDensity, 1.0, 0, 9999.0)); 
  
  //no non-dimensionalization default
  paramListReal.push_back(Parameter<Type>("refPressure", &this->ref_pressure, REF_PRESSURE_DEFAULT, 0.0, 300000.0));
  paramListReal.push_back(Parameter<Type>("refDensity", &this->ref_density, REF_DENSITY_DEFAULT, 0.0, 300.0));
  paramListReal.push_back(Parameter<Type>("refVelocity", &this->ref_velocity, 1.0, 0.0, 9999.0));
  paramListReal.push_back(Parameter<Type>("refLength", &this->ref_length, 1.0, 0.0, 9999999.0));
  paramListReal.push_back(Parameter<Type>("refThermalConductivity", &this->ref_k, 1.0, 0.0, 99999.0));
  paramListReal.push_back(Parameter<Type>("refViscosity", &this->ref_viscosity, 1.0, 0.0, 99999.0));
  paramListReal.push_back(Parameter<Type>("gaussianX", &this->gaussianXloc, 0.0, -999.0, 999.0));
  paramListReal.push_back(Parameter<Type>("gaussianY", &this->gaussianYloc, 0.0, -999.0, 999.0));
  paramListReal.push_back(Parameter<Type>("gaussianAmplitude", &this->gaussianAmpl, 0.25, -999.0, 999.0));
  paramListReal.push_back(Parameter<Type>("gaussianVelocityAmplitude", &this->gaussianVelAmpl, 0.00, 0.0, 999.0));
  paramListReal.push_back(Parameter<Type>("gravityMagnitude", &this->gravity, 9.80665, 0.0, 999.0));

  //DESIGN SECTION
  this->mode = 0;
  paramListInt.push_back(Parameter<Int>("designSolver", &this->designSolver, 0, 0, 1)); //default to sgs
  paramListInt.push_back(Parameter<Int>("designSGSIterations", &this->designNsgs, 1500, 0, 9999));
  paramListInt.push_back(Parameter<Int>("preconditioner", &this->designPrecond, 0, 0, 10)); //no preconditioner
  paramListInt.push_back(Parameter<Int>("designRestarts", &this->designRestarts, 1, 1, 9999)); //default no restarted GMRES
  paramListInt.push_back(Parameter<Int>("designSearchDirections", &this->designSearchDir, 0, 0, 99999)); //default no GMERS
  paramListInt.push_back(Parameter<Int>("objectiveFunction", &this->objFuncId, 0, 0, 99999));
  paramListInt.push_back(Parameter<Int>("sensorEqn", &this->sensEqn, 0, 0, 9999));

}

template <class Type>
Int Param<Type>::ReadSpace(std::string spaceName, std::ifstream& fin, 
			   std::streampos locBegin, std::streampos locEnd)
{
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;

  spacename = spaceName;
  casestring = path + spaceName;
  fin.clear();
  fin.seekg(locBegin);
  std::cout << "PARAM: Reading solution space \"" << spaceName << "\" from parameter file" << std::endl;
  while(fin.tellg() < locEnd){
    c = fin.peek();
    if(c == '#' || c == ' ' || c == '\n'){
      getline(fin, trash);
      trash.clear();
      continue;
    }
    else{
      getline(fin, temp);
      err = ParseLine(temp);
      if(err != 0){
	Abort << "PARAM: could not parse line -- " + temp +  
	  "\nPARAM: check that correct keyword is used";
      }
      err2 += err;
      temp.clear();
    }
  }
  return (err2);
}

template <class Type>
Int Param<Type>::ParseLine(std::string& line)
{
  if(line.length() == 0){
    return (0);
  }

  //Attempt to parse the line with every keyword until we get a successful
  //parse or try them all
  for(std::vector<ParameterString>::iterator it = paramListString.begin();
      it != paramListString.end(); ++it){
    ParameterString& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(std::vector<ParameterStringList>::iterator it = paramListStringList.begin();
      it != paramListStringList.end(); ++it){
    ParameterStringList& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(std::vector<ParameterBool>::iterator it = paramListBool.begin();
      it != paramListBool.end(); ++it){
    ParameterBool& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(std::vector<Parameter<Int> >::iterator it = paramListInt.begin();
      it != paramListInt.end(); ++it){
    Parameter<Int>& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(typename std::vector<Parameter<Type> >::iterator it = paramListReal.begin();
      it != paramListReal.end(); ++it){
    Parameter<Type>& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(typename std::vector<ParameterList<Type> >::iterator it = paramListRealList.begin();
      it != paramListRealList.end(); ++it){
    ParameterList<Type>& p = *it;
    if(p.ParseLine(line)) return 0;
  }
  for(std::vector<ParameterEnum >::iterator it = paramListEnum.begin();
      it != paramListEnum.end(); ++it){
    ParameterEnum& p = *it;
    if(p.ParseLine(line)) return 0;
  }


  return (1);
}

template <class Type>
void Param<Type>::PrintAllParams()
{
 for(std::vector<ParameterString>::iterator it = paramListString.begin();
      it != paramListString.end(); ++it){
    ParameterString& p = *it;
    p.Print();
  }
  for(std::vector<ParameterStringList>::iterator it = paramListStringList.begin();
      it != paramListStringList.end(); ++it){
    ParameterStringList& p = *it;
    p.Print();
  }
  for(std::vector<ParameterBool>::iterator it = paramListBool.begin();
      it != paramListBool.end(); ++it){
    ParameterBool& p = *it;
    p.Print();
  }
  for(std::vector<Parameter<Int> >::iterator it = paramListInt.begin();
      it != paramListInt.end(); ++it){
    Parameter<Int>& p = *it;
    p.Print();
  }
  for(typename std::vector<Parameter<Type> >::iterator it = paramListReal.begin();
      it != paramListReal.end(); ++it){
    Parameter<Type>& p = *it;
    p.Print();
  }
  for(typename std::vector<ParameterList<Type> >::iterator it = paramListRealList.begin();
      it != paramListRealList.end(); ++it){
    ParameterList<Type>& p = *it;
    p.Print();
  }
  for(std::vector<ParameterEnum >::iterator it = paramListEnum.begin();
      it != paramListEnum.end(); ++it){
    ParameterEnum& p = *it;
    p.Print();
  }
}

template <class Type>
void Param<Type>::PostCompute()
{
  ref_time = ref_length/ref_velocity;
  //Type v3 = ref_velocity*ref_velocity*ref_velocity;
  //ref_k = ref_density*v3*ref_length/ref_temperature;
  //ref_viscosity = ref_density*ref_velocity*ref_length;
  ref_enthalpy = ref_velocity*ref_velocity;

  // warning: only one should ever be specified at a time
  // we assume reference pressure takes precedence if both do show up in the input deck
  if(real(ref_pressure) != REF_PRESSURE_DEFAULT){
    ref_density = ref_pressure/(ref_velocity*ref_velocity);
  }
  else if(real(ref_density) != REF_DENSITY_DEFAULT){
    ref_pressure = ref_density*ref_velocity*ref_velocity;
  }

  Type Rs = R_UNIV/(MW/1000.0); //J/kg.K
  std::cout << "PARAM: specific gas constant - " << Rs << " (J/kg.K)" << std::endl;
  
  Type Cv = Cp - Rs; //Cp = Cv + Rs - for ideal gas only
  std::cout << "PARAM: specific heat Cv - " << Rs << " (J/kg.K)" << std::endl;
  this->gamma = Cp/Cv;

  //Compute the required reference temperature to ensure that the isentropic relations hold
  // i.e. (T*) = gamma * (P*) / (rho*) and that the ideal gas equations hold dimensionally T = P/(rho*Rs)
  // and (T*) * Tref = T
  // this implies  Tref = {P/(rho*Rs)} / {gamma * (P*) /(rho*)} = Pref/(gamma * rhoref * Rs)
  ref_temperature = ref_pressure/(gamma * ref_density * Rs);
  
  //Reynolds number is now set in eqnset object initialization
  if(eqnset_id == IncompressibleNS || eqnset_id == CompressibleNS || 
     eqnset_id == CompressibleNSFR){
    viscous = 1;
  }
  
  if(limiterRefresh <= 0){
    limiterRefresh = 1;
  }

  if(sorder > 1){
    if(sorder == 2){
      chi = 0.0;      
    }
    if(sorder >= 3){
      chi = 0.5;
    }
  }
  
  if(fieldJacType == 2 || boundaryJacType == 2){
    requireComplex = true;
  }
  if(dynamicCFL){
    cflPrev = cflStart;
  }
  //if the chemical database path has been set add the path
  //we assume the user means a local file
  if(chemDB != defaultChemDB){
    chemDB = path + chemDB;
  }

  //make sure both molefractions and massfractions are not set
  if(molefractions.size() != 0 && massfractions.size() != 0){
    Abort << "WARNING: mass fractions and mole fraction cannot both be set in input file";
  }

}

template <class Type>
void Param<Type>::PrintSolverParams()
{

  std::cout << "/******************************************/\n";
  std::cout << "  SOLUTION SPACE: " << spacename << std::endl;
  std::cout << "/******************************************/\n\n";
    
  std::cout << "PARAM: Runtime parameters from file" << std::endl;
  std::cout << "===================================" << std::endl;
  std::cout << "\tEquation set: " << eqnSets[eqnset_id] << std::endl;
  std::cout << "\tFlux type: " << fluxTypes[flux_id] << std::endl;
  if(eqnset_id == IncompressibleEuler || eqnset_id == IncompressibleNS){
    std::cout <<"\tArtificial compressibility (beta): " << beta << std::endl;
  }
  if(viscous){
    std::cout << "\tViscosity terms are on!" << std::endl;
    std::cout << "\t\tPrandtl number: " << Pr << std::endl;
    std::cout << "\t\tTurbulent Prandtl number: " << PrT << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tTurbulence model: " << turbModel << std::endl;
    std::cout << "\t\tTurbulence model spatial order: " << turbModelSorder << std::endl;
  }
  if(gravity_on){
    std::cout << "\tGravity buoyancy terms ON" << std::endl;
    std::cout << "\t\tGravity vector: " << gravdir[0] << " " << gravdir[1] << " " << gravdir[2] << std::endl;
    std::cout << "\t\tGravity magnitude (dimensional): " << gravity << std::endl;
  }
  else{
    std::cout << "\tGravity buoyancy terms OFF" << std::endl;
  }
  if(nSgs <= 0){
    std::cout << "\tUsing explicit solve -- nSgs <= 0" << std::endl;
  }
  else{
    std::cout << "\tNumber of SGS sweeps: " << nSgs << std::endl;
    std::cout << "\tJacobian update frequency: " << jacobianFreq << std::endl;
    std::cout << "\tField jacobian type: ";
    switch(fieldJacType){
    case 0:
      std::cout << "upwind" << std::endl;
      break;
    case 1:
      std::cout << "central difference" << std::endl;
      break;
    case 2:
      std::cout << "complex" << std::endl;
      break;
    default:
      std::cout << "upwind" << std::endl;
      break;
    }
    std::cout << "\tBoundary jacobian type: ";
    switch(boundaryJacType){
    case 0:
      std::cout << "upwind" << std::endl;
      break;
    case 1:
      std::cout << "central difference" << std::endl;
      break;
    case 2:
      std::cout << "complex" << std::endl;
      break;
    default:
      std::cout << "upwind" << std::endl;
      break;
    }
    if(boundaryJacEval){
      std::cout << "\tWARNING: Using approximate boundary jacobians!!!" << std::endl;
    }
    if(hojac){
      std::cout << "\tUsing higher order jacobians" << std::endl;
    }
  }
  std::cout << "\tTemporal order of accuracy: " << torder << std::endl;
  std::cout << "\tSpatial order of accuracy: " << sorder << std::endl;
  if(sorder > 1){
    std::cout << "\tNumber of first order steps to take: " << nFirstOrderSteps << std::endl;
    if(limiter <= 0){
      std::cout << "\tLimiter OFF" << std::endl;
    }
    else{
      std::cout << "\tLimiter ON" << std::endl;
      std::cout << "\t\tLimter refresh rate: " << limiterRefresh << std::endl;
      switch(limiter){
      case 1:
	std::cout << "\t\tUsing Barth limiter" << std::endl;
	break;
      case 2:
	std::cout << "\t\tUsing Venkatakrishnan limiter" << std::endl;
	break;
      case 3:
	std::cout << "\t\tUsing modified Venkatakrishnan limiter" << std::endl;
	break;
      default:
	std::cout << "\t\tWARNING: Limiter not found" << std::endl;
	std::cout << "\t\tDefaulting to Barth!!!" << std::endl;
	break;
      }
    }
  }
  if(useLocalTimeStepping){
    std::cout << "\tUsing Local Timestepping!!" << std::endl;
    std::cout << "\tCFL: " << cfl << std::endl;
    std::cout << "\tNondimensional timestep size: " << dt << std::endl;
    std::cout << "\tDimensional timestep size: " << dt*ref_time << std::endl;
    //set previous iteration cfl equal to cfl
    if(cflRampingSteps > 0){
      std::cout << "\tCFL ramping enabled!!" << std::endl;
      std::cout << "\tCFL start: " << cflStart << std::endl;
      std::cout << "\tCFL ramping steps: " << cflRampingSteps << std::endl;
    }
    if(dynamicCFL){
      std::cout << "\tDynamic CFL enabled!!" << std::endl;
    }
    if(real(dt) > 0.0){
      std::cout << "\tRunning unsteady" << std::endl;
    }
    else{
      std::cout << "\tRunning steady" << std::endl;
    }
  }
  else{
    std::cout << "\tUsing Minimum Timestepping!!" << std::endl;
    if(real(dt) > 0.0){
      std::cout << "\tRunning unsteady" << std::endl;
      std::cout << "\tNondimensional timestep size: " << dt << std::endl;
      std::cout << "\tDimensional timestep size: " << dt*ref_time << std::endl;
    }
    else{
      std::cout << "Running steady" << std::endl;
      std::cout << "\tCFL: " << cfl << std::endl;
      if(cflRampingSteps > 0){
	std::cout << "\tCFL ramping enabled!!" << std::endl;
	std::cout << "\tCFL start: " << cflStart << std::endl;
	std::cout << "\tCFL ramping steps: " << cflRampingSteps << std::endl;
      }
    }
  }
  std::cout << "\tVelocity: " << velocity << std::endl;
  if(velocityRampingSteps > 0){
    std::cout << "\tVelocity ramping enabled!!" << std::endl;
    std::cout << "\tVelocity start: " << velocityStart << std::endl;
    std::cout << "\tVelocity ramping steps: " << velocityRampingSteps << std::endl;
  }
  if(rxnOn){
    std::cout << "Chemical reactions turned ON" << std::endl;
  }
  else{
    std::cout << "Chemical reactions turned OFF" << std::endl;
  }
  if(no_cvbc){
    std::cout << "CVBCs turned OFF" << std::endl;
  }
  std::cout << "Reference variables (specified OR derived)" << std::endl;
  std::cout << "===============================" << std::endl;
  std::cout << "\tvelocity : " << ref_velocity << std::endl;
  std::cout << "\tlength : " << ref_length << std::endl;
  std::cout << "\ttemperature : " << ref_temperature << std::endl;
  std::cout << "\tpressure : " << ref_pressure<< std::endl;
  std::cout << "\tviscosity: " << ref_viscosity << std::endl;
  std::cout << "\tconductivity: " << ref_k << std::endl;
  std::cout << "\ttime : " << ref_time << std::endl;
  std::cout << "\tenthalpy: " << ref_enthalpy << std::endl;
  std::cout << "\tdensity: " << ref_density << std::endl;
  std::cout << "\tGamma (Cp/Cv): " << gamma << std::endl;
  std::cout << "\tMW: " << MW << " (g/mol)" << std::endl;
  std::cout << "\tCp: " << Cp << std::endl;
  std::cout << std::endl;
  std::cout << "\tFlow direction vector " << flowdir[0] <<" "<< flowdir[1] <<" "<< flowdir[2] << std::endl;
  std::cout << "\tLift direction vector " << liftdir[0] <<" "<< liftdir[1] <<" "<< liftdir[2] << std::endl;
  std::cout << std::endl;
  if(useRestart){
    std::cout << "\tAttempting to use restarted flow solution" << std::endl;
    if(preserveRestartCounters){
      std::cout << "\t\tPreserving iteration counters" << std::endl;
    }
    else{
      std::cout << "\t\tNot preserving iteration counters" << std::endl;
    }
  }
  if(solutionWrite){
    std::cout << "\tUsing unsteady solution writing every " << solutionWrite << " timesteps" << std::endl;
  }
  if(movement){
    std::cout << "\tUsing moving mesh solution routines" << std::endl;
    if(gcl){
      std::cout << "\t\tGCL enabled" << std::endl;
    }
    else{
      std::cout << "\t\tGCL NOT enabled - WARNING: if mesh deforms, you are in trouble" << std::endl;
    }
  }
  std::cout << std::endl;
  if(mode != 0){
    std::cout << "Design parameters" << std::endl;
    std::cout << "=================" << std::endl;
    std::cout << "Design mode: ";
    if(mode == 1){
      std::cout << "Objective function evaluation";
    }
    else if(mode == 2){
      std::cout << "Direct mode design sensitivity";
    }
    else if(mode == 3){
      std::cout << "Adjoint mode design sensitivity";
    }
    else if(mode == 4){
      std::cout << "CTSE mode design sensitivity";
    }
    else if(mode == 5){
      std::cout << "Mesh smoothing";
    }
    else if(mode == 6){
      std::cout << "Mesh sensitivity evaluation";
    }
    std::cout << std::endl;
    if(designSolver == 0){
      std::cout << "Solver type: Using SGS solver" << std::endl;
      std::cout << "Number of iterations: " << designNsgs << std::endl;
    }
    else if(designSolver == 1){
      std::cout << "Solver type: Using GMRES solver" << std::endl;
      std::cout << "Number of search directions: " << designSearchDir << std::endl;
      std::cout << "Number of restarts: " << designRestarts << std::endl;
    }
    std::cout << "Preconditioner: " << designPrecond << std::endl; 
  }
  if(gaussianSource){
    std::cout << "Applying gaussian source term to equation " << gaussianEqn << std::endl;
    std::cout << "Applying gaussian source to BC " << gaussianBCid << std::endl;
    std::cout << "Gaussian source amplitude " << gaussianAmpl << std::endl;
    std::cout << "Gaussian velocity amplitude " << gaussianVelAmpl << std::endl;
    std::cout << "Gaussian x location " << gaussianXloc << std::endl;
    std::cout << "Gaussian y location " << gaussianYloc << std::endl;
  }
  std::cout << std::endl;

  return;
}

template <class Type>
Type Param<Type>::GetCFL()
{
  return(cflPrev);
}

template <class Type>
void Param<Type>::UpdateCFL(Int iter, Type deltaResidual)
{
  if(dynamicCFL){ //CFL solution steering 
    if(iter > 1){
      //allow 10% change per iteration
      Type delta = 0.1;
      Type sgn;
      if(real(deltaResidual) > 0.0){
	sgn = 1.0;
      }
      else{
	sgn = -1.0;
      }
      //A = MAX(0, previousCFL * delta * sign (+/-)) --- only take positive CFLS
      Type A = MAX(0.0, cflPrev + delta * sgn * cflPrev);
      Type cflnew = MIN(A, cfl);
      cflPrev = cflnew;
    }
    else{
      cflPrev = cflStart;
    }
  }
  else{ //prescribed CFL ramping over time
    if(cflRampingSteps > 0 && iter < cflRampingSteps){
      cflPrev = ( (cfl-cflStart) / (Type)cflRampingSteps * (Type)iter + cflStart);
    }
    else{
      cflPrev = cfl;
    }
  }
}

//This function returns the non-dimensional velocity as non-dimensionalized by ref_velocity
//if your function needs velocity non-dimensionalized by some other value like the speed of
//sound as in compressible equations, that needs to be performed before this value is used
template <class Type>
Type Param<Type>::GetVelocity(Int iter)
{
  if(velocityRampingSteps > 0 && iter < velocityRampingSteps){
    return ( (velocity-velocityStart) / (Type)velocityRampingSteps * (Type)iter + velocityStart);
  }
  else{
    return (velocity);
  }
}

template <class Type>
void Param<Type>::UpdateForDesign()
{
  std::string designName = path+spacename;

  Int ndv = GetNdvDesignFile(designName);
  Type* x = new Type[ndv];
  Type* bounds = new Type[ndv*2];
  Type f;
  Type* grad = new Type[ndv];
  Int* dvType = new Int[ndv];
  Int parType;

  if(ReadDesignFile(designName, &ndv, x, bounds, &f, grad, dvType)){
    Abort << "Could not read design file";
  }

  for(Int beta = 0; beta < ndv; beta++){
    Int parType = dvType[beta];
    if(parType == 0 || parType == 1 || parType == 2){
      //nothing - direct node movement
    }
    else if(parType == 3 || parType == 4 || parType == 5){
      //nothing - hicks-henne function
    }
    else if(parType == 6 || parType == 7 || parType == 8){
      //nothing - movement of whole boundary (x,y,z)
    }
    else if(parType == 9){
      //x location of gaussian plume source
      gaussianXloc = x[beta];
    }
    else if(parType == 10){
      //y location of gaussian plume source
      gaussianYloc = x[beta];
    }
    else{
      //should never be here
      std::stringstream ss;
      ss << "In Param::UpdateForDesign() and design parameter type ";
      ss << parType;
      ss << " not found!";
      Abort << ss.str();
    }
  }

  PostCompute();
  
  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;
}

template <class Type>
const ParameterString& Param<Type>::GetParameterString(std::string name) const
{
  for(std::vector<ParameterString>::const_iterator it = paramListString.begin();
      it != paramListString.end(); ++it){
    const ParameterString& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterString() dead";
  return(*paramListString.begin());
}

template <class Type>
const ParameterStringList& Param<Type>::GetParameterStringList(std::string name) const
{
  for(std::vector<ParameterStringList>::const_iterator it = paramListStringList.begin();
      it != paramListStringList.end(); ++it){
    const ParameterStringList& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterStringList() dead";
  return(*paramListStringList.begin());
}

template <class Type>
const ParameterBool& Param<Type>::GetParameterBool(std::string name) const
{
  for(std::vector<ParameterBool>::const_iterator it = paramListBool.begin();
      it != paramListBool.end(); ++it){
    const ParameterBool& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterBool() dead";
  return(*paramListBool.begin());
}

template <class Type>
const Parameter<Int>& Param<Type>::GetParameterInt(std::string name) const
{
  for(std::vector<Parameter<Int> >::const_iterator it = paramListInt.begin();
      it != paramListInt.end(); ++it){
    const Parameter<Int>& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterInt() dead";
  return(*paramListInt.begin());
}

template <class Type>
const Parameter<Type>& Param<Type>::GetParameterType(std::string name) const
{
  for(typename std::vector<Parameter<Type> >::const_iterator it = paramListReal.begin();
      it != paramListReal.end(); ++it){
    const Parameter<Type>& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterType() dead";
  return(*paramListReal.begin());
}

template <class Type>
const ParameterList<Type>& Param<Type>::GetParameterListType(std::string name) const
{
  for(typename std::vector<ParameterList<Type> >::const_iterator it = paramListRealList.begin();
      it != paramListRealList.end(); ++it){
    const ParameterList<Type>& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterListType() dead";
  return(*paramListRealList.begin());
}

template <class Type>
const ParameterEnum& Param<Type>::GetParameterEnum(std::string name) const
{
  for(std::vector<ParameterEnum >::const_iterator it = paramListEnum.begin();
      it != paramListEnum.end(); ++it){
    const ParameterEnum& p = *it;
    if(p.GetKey() == name){
      return p;
    }
  }
  Abort << "Param::GetParameterEnum() dead";
  return(*paramListEnum.begin());
}

//This is the master param file reading call, will create and load a vector
//of param files which correspond to each solution space defined in the param file
template <class Type> 
Int ReadParamFile(std::vector<Param<Type>*>& paramFileList, std::string casename, std::string pathname, 
		  Bool verbose)
{
  std::ifstream fin;
  std::string trash, temp;
  Int err = 0;
  Int err2 = 0;
  char c;

  std::string endBang = ">>>";
  std::string beginSpace = "<<<BEGIN SPACE ";
  std::string endSpace = "<<<END SPACE>>>";
  std::string filename = pathname+casename + ".param";

  fin.open(filename.c_str());

  if(fin.is_open()){
    std::cout << "PARAM: Reading parameter file --> " << filename << std::endl;
    //parse all of the solution spaces from the file
    while(getline(fin, temp)){
      if(temp.size() != 0){
	c = temp[0];
	if(c == '#' || c == ' ' || c == '\n'){
	  temp.clear();
	  continue;
	}
      }
      else{
	continue;
      }
      size_t loc;
      loc = temp.find(beginSpace);
      if(loc != std::string::npos){
	std::string spaceName;
	if(GetStringBetween(beginSpace, endBang, temp, spaceName)){
	  Abort << "WARNING: Solution space delimiter found uncomplete, given -- " + temp;
	  return(1);
	}
	std::streampos beginSpacepos = fin.tellg();
	std::streampos endSpacepos = fin.tellg();
	//we have the beginning of the solution space defintion, now find the end
	std::streampos preEndSpacepos;
	while(getline(fin, trash)){
	  if(trash.size() != 0){
	    c = trash[0];
	    if(c == '#' || c == ' ' || c == '\n'){
	      trash.clear();
	      continue;
	    }
	  }
	  else{
	    continue;
	  }
	  //this occurs if there is not a newline at the end of the file
	  loc = trash.find(endSpace);
	  if(loc != std::string::npos){
	    endSpacepos = preEndSpacepos + (std::streampos)endSpace.size();
	    Param<Real>* parampt = new Param<Real>();
	    parampt->path = pathname;
	    paramFileList.push_back(parampt);
	    if(paramFileList.back()->ReadSpace(spaceName, fin, beginSpacepos, preEndSpacepos)){
	      return (1);
	    }
	    //seek back to where we were to continue, routine readspace may modify stream
	    fin.seekg(endSpacepos);
	    break;
	  }
	  //get the position before the next line is read, might be the end of space delimiter
	  preEndSpacepos = fin.tellg();
	}
	if(beginSpacepos == endSpacepos){
	  Abort << "WARNING: Solution space end delimiter not found";
	  return(1);
	}
	trash.clear();
      }
      temp.clear();
    }
    fin.close();
  }
  else{
    Abort << "PARAM READFILE: Cannot open param file --> " + filename;
    return (1);
  }

  if(paramFileList.size() == 0){
    Abort << "PARAM READFILE: Did not find any defined solution spaces\n\tThese are defined by <<<BEGIN SPACE (name)>>> and <<<END SPACE>>> delimiters";
    return (1);
  }
  

  for(typename std::vector<Param<Real>*>::iterator it = paramFileList.begin(); it != paramFileList.end(); ++it){
    (*it)->PostCompute();
  }

  return (err2);
}
