#ifndef VARS_ID_H__
#define VARS_ID_H__

//WARNING: make sure you know what you are doing before changing the ordering 
//here... dependencies scattered around

//enumeration for computed variables
enum VarsTypes
  {
    WallDistance,
    LimiterValue,
    TurbulentViscosity,
    SurfaceCp,
    SurfaceYp,
    SurfaceCf,
    SurfaceQdot,
  };

const std::string VarsNames[] = 
  {
    "wallDistance",
    "limiterValue",
    "turbulentViscosity",
    "surfaceCp",
    "surfaceYp",
    "surfaceCf",
    "surfaceQdot",
  };


#endif
