#ifndef FFD_H__
#define FFD_H__

#include "general.h"
#include "mesh.h"
#include "bc.h"
#include <vector>

//forward declarations
template <class Type> class BoundaryConditions;
template <class Type> class Mesh;

/*
This class is design to perform FFD (free form deformation) of a boundary which
is identified through input.  The end goal here being the ability to deform a mesh
in a smooth and continuously differentiable way for design optimization purposes.
*/
template <class Type>
class FFD
{
public:
  template <class Type2>
  FFD(const Int boundaryAttached, const Mesh<Type>* m, const BoundaryConditions<Type2>* bcs);
  ~FFD();

protected:

private:
  FFD();

  Int bcAttached; //Boundary id to which FFD box is attached
  Type xmax;  //max x-extent
  Type xmin;  //min x-extent
  Type ymax;  //max y-extent
  Type ymin;  //min y-extent
  Type zmax;  //max z-extent
  Type zmin;  //min z-extent

};

#endif
