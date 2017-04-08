#include "etypes.h"

std::string ETypeToName(Int type)
{
  switch(type){
  case TRI:
    return "Tri";
  case QUAD:
    return "Quad";
  case TET:
    return "Tet";
  case PYRAMID:
    return "Pyramid";
  case PRISM:
    return "Prism";
  case HEX:
    return "Hex";
  default:
    return "Bad Element Type";
  }
};
