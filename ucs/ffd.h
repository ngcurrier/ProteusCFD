#ifndef FFD_H__
#define FFD_H__

#include "general.h"
#include "mesh.h"
#include "bc.h"
#include <vector>

/*
This class is design to perform FFD (free form deformation) of a boundary which
is identified through input.  The end goal here being the ability to deform a mesh
in a smooth and continuously differentiable way for design optimization purposes.
*/
template <class Type>
class FFD
{
public:
  FFD(const int boundaryAttached, const Mesh* m, const BoundaryConditions* bcs);
  ~FFD();

protected:

private:
  FFD();

  int bcAttached; //Boundary id to which FFD box is attached
  double xmax;  //max x-extent
  double xmin;  //min x-extent
  double ymax;  //max y-extent
  double ymin;  //min y-extent
  double zmax;  //max z-extent
  double zmin;  //min z-extent

};

#endif
