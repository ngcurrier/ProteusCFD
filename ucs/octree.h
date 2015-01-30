#ifndef OCTREE_H__
#define OCTREE_H__

#include "general.h"
#include "macros.h"
#include <vector>
#include <set>

template <class Type> class BBox;

template <class Type>
class Octree
{
public:
  Octree(Type* xyz, Int npts, Int maxDepth, Int minNodesPerLevel);
  ~Octree();

  //given a new point xyz, returns the distance to the nearest node
  Type FindDistance(Type* xyz, Int& nodeId);

  //contracts each bounding box to its minimum extents, makes searching more accurate
  void Contract();

  //sets the six coordinates (xmin, xmax, ymin, ymax, zmin, zmax) for the
  //root bounding box
  void GetExtents(Type* box);

  void Print(std::ostream& str);

  std::vector<BBox<Type>*> flatList;
  BBox<Type>* root;
  Type* xyz;
  Int deepestLevel;
private:
  Int maxdepth;
  Int minNodes;
  Int npts;
};

//class bounding box, used by octree class to find points in a field
template <class Type>
class BBox
{
public:
  BBox();
  ~BBox();

  //level id
  Int level;

  //extents
  Type max[3];
  Type min[3];

  //number of points within this level
  Int count;

  //number of boxes which live below this one, we prune zero
  //occupancy boxes b/c of storage requirements for large meshes
  Int nBelow;

  //list of points within this level
  //note this will reference original list 
  std::vector<Int> list;

  //for octree structure, pointer down a level
  //since each split contains eight children
  //this is a list of length 8
  BBox<Type>** Down;

  //store a reference to the octree the box is a part of
  Octree<Type>* octree;

  //add a level to the octree
  void Split();
  //check to see if this node is split
  Bool IsSplit();
  //build the list of nodes contained within self
  void BuildList(BBox<Type>* Up);
  //returns the distance to the nearest node in the octree, and the id in the xyz list
  Type GetNodeDistance(Type* xyz, Int& nodeId);
  //returns pointer to voxel containing that node, also passes back distance to leaf
  BBox<Type>* GetLeaf(Type* xyz, Type* dist);
  //returns the edge distance to the voxel
  Type GetVoxelDistance(const Type* xyz);
  //returns the distance from voxel centroid 
  Type GetVoxelDistanceCentroid(const Type* xyz);
  //Returns the six coordinates (xmin, xmax, ymin, ymax, zmin, zmax) for the bounding box
  void GetExtents(Type* box);
  //Contracts each bounding box to the minimum extents of the nodes it contains
  void Contract();
  //returns a vector with pointers to the bounding boxes which are neighbors of a leaf
  //that is closer than our current best guess
  void GetNeighbors(std::vector<BBox<Type>*>& neighbors, const Type* xyz, const Type dist);
private:

};

//include implementation
#include "octree.tcc"

#endif
