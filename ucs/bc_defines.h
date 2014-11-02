#ifndef BC_DEFINES_H__
#define BC_DEFINES_H__

enum BCTypes
  {
    //do not move ParallelBoundary from first position
    //or the house of cards in bc.cpp will come crashing down
    ParallelBoundary,
    Dirichlet,
    Neumann,
    ImpermeableWall,
    NoSlip,
    FarFieldViscous,
    FarField,
    SonicInflow,
    SonicOutflow,
    Symmetry,
    InternalInflowDuctBL,
    InternalInflow,
    InternalOutflow,
    NormalInFarField,
    NormalOutFarField,
    PitchingFarField,
    TotalTempAndPressure,
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
    "totalTempAndPressure"
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
    "twall"
  };

#endif
