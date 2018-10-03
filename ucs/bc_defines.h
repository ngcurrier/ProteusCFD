#ifndef BC_DEFINES_H__
#define BC_DEFINES_H__

enum BCTypes
  {
    //do not move ParallelBoundary from first position
    //or the house of cards in bc.cpp will come crashing down
    Proteus_ParallelBoundary,
    Proteus_Dirichlet,
    Proteus_Neumann,
    Proteus_ImpermeableWall,
    Proteus_NoSlip,
    Proteus_FarFieldViscous,
    Proteus_FarField,
    Proteus_SonicInflow,
    Proteus_SonicOutflow,
    Proteus_Symmetry,
    Proteus_InternalInflowDuctBL,
    Proteus_InternalInflow,
    Proteus_InternalOutflow,
    Proteus_NormalInFarField,
    Proteus_NormalOutFarField,
    Proteus_PitchingFarField,
    Proteus_TotalTempAndPressure,
    Proteus_HeatFlux,
    Proteus_Isothermal,
    Proteus_PythonBC,
    Proteus_NULL,
    NUM_BC_TYPES
  };
const std::string BCs[] =
  {
    "parallelBoundary-DoNotUse-DoNotChange",
    "dirichlet",
    "neumann",
    "impermeableWall",
    "noSlip",
    "farFieldViscous",
    "farField",
    "sonicInflow",
    "sonicOutflow",
    "symmetry",
    "internalInflowDuctBL",
    "internalInflow",
    "internalOutflow",
    "normalInFarField", //normal inflow to flat boundary, subject to farfield conditions
    "normalOutFarField",  //normal ouflow to flat boundary, subject to farfield conditions
    "pitchingFarField",
    "totalTempAndPressure",
    "heatFlux",
    "isothermal",
    "pythonBC",
    "NULL"
  }; 

enum BCvarsId
  {
    QREF,
    FlowDirection,
    Velocity,
    BackPressure,
    Density,
    Axis,
    Point, 
    Omega,
    SlipSpeed,
    SlipDirection,
    MassFractions,
    Moving,
    BleedSteps,
    Twall,
    Flux,
    NUM_BC_VAR_IDS
  };

const std::string BCVars[] =
  {
    "Qref",
    "flowDirection",
    "velocity",
    "backPressure",
    "density",
    "axis",
    "point",
    "omega",
    "slipSpeed",
    "slipDirection",
    "massFractions",
    "isMoving",
    "bleedSteps",
    "twall",
    "flux"
  };

#endif
