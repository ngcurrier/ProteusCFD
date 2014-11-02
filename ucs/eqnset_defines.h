#ifndef EQNSET_DEFINES_H__
#define EQNSET_DEFINES_H__

#include <string>

////////////////////////////////////////////////////////
//Definitions for param file reader... change
//with high caution -- likely to break 
///////////////////////////////////////////////////////
enum EqnTypes
  {
    CompressibleEulerFR,
    CompressibleNSFR,
    CompressibleEuler,
    CompressibleNS,
    IncompressibleEuler,
    IncompressibleNS,
    HeatTransfer,
    NUM_EQN_SETS
  };
const std::string eqnSets[] = 
  {
    "compressibleEulerFR",
    "compressibleNSFR",
    "compressibleEuler",
    "compressibleNS",
    "incompressibleEuler",
    "incompressibleNS",
    "heatTransfer",
    //this matches EqnTypes enumeration in eqnset.h
    //must be consistent --- be warned!!!
  };

//Type identifier added to avoid conflict with member function name
enum FluxTypes
  {
    AlgebraicFluxType,
    RoeFluxType,
    HLLCFluxType,
    NUM_FLUX_TYPES
  };
const std::string fluxTypes[] = 
  {
    "algebraicFlux",
    "roeFlux",
    "HLLCFlux"
    //this matches FluxTypes enumeration in eqnset.h
    //must be consistent --- be warned!!!
  };
#endif
