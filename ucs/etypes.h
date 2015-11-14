#ifndef ETYPE_H__
#define ETYPE_H__

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

#endif
