#ifndef ETYPE_H__
#define ETYPE_H__

#include <string>
#include "general.h"

//define the element types for the .su2 mesh format
#define SU2_LINE 3
#define SU2_TRI 5
#define SU2_QUAD 9
#define SU2_TET 10
#define SU2_HEX 12
#define SU2_PRISM 13
#define SU2_PYRAMID 14

//Element type hierarchy: these are primitives 
//                        read from mesh file only
#define TRI 0
#define QUAD 1
#define TET 2
#define PYRAMID 3
#define PRISM 4
#define HEX 5
#define MAX_E_TYPES 6

std::string ETypeToName(Int type);

//define the element types for the gmsh mesh format
#define GMSH_LINE 1 
#define GMSH_TRI 2
#define GMSH_QUAD 3
#define GMSH_TET 4
#define GMSH_HEX 5
#define GMSH_PRISM 6
#define GMSH_PYRAMID 7
#define GMSH_POINT 15

#endif
