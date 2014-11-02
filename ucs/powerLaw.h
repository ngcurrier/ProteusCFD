#ifndef POWER_LAW_H__
#define POWER_LAW_H__

#include "general.h"
#include "macros.h"
#include <iostream>
#include <cmath>

//assumes a distance from the leading edge of a flatplate
//takes a walldistance and uinf and returns the powerlaw value of streamwise velocity u
template <class Type, class Type2>
Type PowerLawU(Type uinf, Type wallDist, Type2 Re)
{
  Type re = Re;
  //distance from leading edge
  Type x = 10.0;
  //thickness of the boundary layer, theoretical results
  Type deltaTurb = 0.382*x/(pow(re, 0.2));
  Type deltaLam = 4.91*x/(sqrt(re));

  Type delta = deltaTurb;

  Type u = uinf*pow(wallDist/delta, 1.0/7.0);
  if(real(u) < real(uinf)){
    return u;
  }
  return(uinf);
}

#endif
