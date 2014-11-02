#ifndef UNSBASE_H__
#define UNSBASE_H__

#include "general.h"

//
// Edge data structure: this structure is used to store
//                      data relevant to the solution procedure, not 
//                      geometry
//
//                      edge always point from low to high vertex index
template <class Type>
class Edges
{
 public:
  Int n[2];              //Two vertices which define the edge
  Type a[4];             //Area vector for the edge a[3] is magnitude
                         //all other components are normalized
  Type wl, wr;           //edge weights for least squares gradient calculation
 private:
};

//
// Half Edge (boundary) data structure: used to store relevant information
//                                      concerning a portion of a boundary
//                                      element
template <class Type>
class HalfEdges
{
 public:
  Int n[2];             //The vertex which this half edge is associated with
  Type a[4];            //Area vector for the half edge a[3] is magnitude
                        //all other components are stored normalized
  Int elem;             //element which half edge is derived from
  Int factag;           //factag acquired from root element
 private:
};


#endif
