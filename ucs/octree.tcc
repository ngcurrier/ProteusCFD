#include <deque>

template <class Type>
Octree<Type>::Octree(Type* xyz, Int npts, Int maxDepth, Int minNodesPerLevel) :
  deepestLevel(0), maxdepth(maxDepth), minNodes(minNodesPerLevel), npts(npts), xyz(xyz)
{
  root = new BBox<Type>;
  root->octree = this;
  //the root level contains all of the points
  root->count = npts;
  flatList.reserve(maxDepth/2 *8);
  //put the root octree voxel in the list
  flatList.push_back(root);
  //set the list of nodes equal to the list of all the nodes passed in
  root->list.reserve(root->count);
  for(Int i = 0; i < npts; i++){
    root->list.push_back(i);
  }
  //find max and min of root octree level
  root->max[0] = root->max[1] = root->max[2] = -1.0e99;
  root->min[0] = root->min[1] = root->min[2] = 1.0e99;
  
  for(Int i = 0; i < npts; i++){
    root->max[0] = MAX(root->max[0], xyz[i*3 + 0]);
    root->max[1] = MAX(root->max[1], xyz[i*3 + 1]);
    root->max[2] = MAX(root->max[2], xyz[i*3 + 2]);
    root->min[0] = MIN(root->min[0], xyz[i*3 + 0]);
    root->min[1] = MIN(root->min[1], xyz[i*3 + 1]);
    root->min[2] = MIN(root->min[2], xyz[i*3 + 2]);
  }

  //force squareness of the root node
  Type maxdim = 0.0;
  for(Int i = 0; i < 3; i++){
    maxdim = MAX(maxdim, root->max[i] - root->min[i]);
  }
  for(Int i = 0; i < 3; i++){
    root->max[i] = root->min[i] + maxdim;
  }
  

  //we don't do tail recursion for fear of blowing the stack...
  BBox<Type>* active = root;
  Int n = 0;
  do{
    if(active->list.size() > minNodes && active->level < maxDepth){
      active->Split();
    }
    //this count and list gets updated inside the split() routine
    //since all new voxels are appended we just keep going until
    //we are out of things to check for splitting conditions
    if(n < flatList.size()){
      active = flatList[n];
    }
    else break;
    n++;
  }while (n < flatList.size());
  
  //contract each bounding box to its minimum size, makes searching easier
  Contract();

  return;
}

template <class Type>
void Octree<Type>::Print(std::ostream& str)
{
  Int mbyte = 1048576;  //since we are referring to memory
  str << "OCTREE: Number of bounding boxes created " << flatList.size() << std::endl;
  str << "OCTREE: Bounding box size: " << sizeof(BBox<Type>) << " bytes" << std::endl;
  str << "OCTREE: Deepest level " << deepestLevel << std::endl;
  str << "OCTREE: size ~= " << (double)(sizeof(BBox<Type>)*flatList.size() + 
					sizeof(BBox<Type>*)*flatList.capacity() + 
					sizeof(Int)*root->count)/(double)mbyte
      << " MB for storage" << std::endl;

  Int mostNodes = 0;
  for(Int i = 0; i < flatList.size(); i++){
    BBox<Type>* box = flatList[i];
    //only check the terminating leafs
    if(!box->IsSplit()){
      if(box->list.size() > mostNodes){
	mostNodes = box->list.size();
      }
    }
  }
  str << "OCTREE: Most nodes in leaves " << mostNodes << std::endl;

  Int empties = 0;
  for(Int i = 0; i < flatList.size(); i++){
    BBox<Type>* box = flatList[i];
    if(!box->IsSplit()){
      if(box->count == 0){
	empties++;
      }
    }
  }
  str << "OCTREE: Number of empty leaves " << empties << std::endl;
}

template <class Type>
Octree<Type>::~Octree()
{
  for(Int i = 0; i < flatList.size(); i++){
    delete flatList[i];
  }
}

template <class Type>
Type Octree<Type>::FindDistance(Type* xyz, Int& nodeId)
{
  return (root->GetNodeDistance(xyz, nodeId));
}

template <class Type>
void Octree<Type>::GetExtents(Type* box)
{
  box[0] = root->min[0];
  box[1] = root->max[0];
  box[2] = root->min[1];
  box[3] = root->max[1];
  box[4] = root->min[2];
  box[5] = root->max[2];
}

template <class Type>
void Octree<Type>::Contract()
{
  std::vector<Bool> contracted(flatList.size(), false);

  Int deepestLevel = 0;

  //contract the deepest nodes first which contain the nodes
  for(Int i = 0; i < flatList.size(); i++){
    BBox<Type>* box = flatList[i];
    if(!box->IsSplit()){
      box->Contract();
      contracted[i] = true;
      if(deepestLevel < box->level) deepestLevel = box->level;
    }
  }
  deepestLevel--;

  //work up starting from the bottom contracting as we go
  while(deepestLevel > 0){
    for(Int i = 0; i < flatList.size(); i++){
      BBox<Type>* box = flatList[i];
      if(box->count != 0 && box->level == deepestLevel && !contracted[i]){
	box->Contract();
	contracted[i] = true;
      }
    }
    deepestLevel--;
  }

}

template <class Type>
void BBox<Type>::GetNeighbors(std::vector<BBox<Type>*>& neighbors, const Type* xyz, const Type dist)
{
  //The game here is to recurse up the tree and add any bounding
  //box that is closer to our point to a queue'ing list
  //since we know the tree isn't balanced, and we'll eventually
  //end up back at the root, we just start at the top and descend

  neighbors.clear();
  
  //we're gonna use a FIFO queue
  std::deque<BBox<Type>*> queue;
  queue.push_back(octree->root);
  std::set<BBox<Type>*> unique;

  do{
    BBox<Type>* current = queue.front();
    queue.pop_front();
    if(current->IsSplit()){
      for(Int i = 0; i < current->nBelow; i++){
	if(current->Down[i]->count != 0){
	  //if this box is closer than the closest point we've found yet
	  //add it to the queue to recurse into
	  if(real(current->Down[i]->GetVoxelDistance(xyz)) < real(dist)){
	    queue.push_back(current->Down[i]);
	  }
	}
      }
    }
    //the current box is a terminating node, if it is closer than the closest
    //point found thus far, it mandates a total search internally
    else{
      if(current->count != 0){
	if(real(current->GetVoxelDistance(xyz)) < real(dist)){
	  unique.insert(current);
	}
      }
    }
  }while(queue.size() != 0);
  
  //copy data over to the neighbors vector since it has random access iterator
  neighbors.resize(unique.size());
  typename std::vector<BBox<Type>*>::iterator vit = neighbors.begin();
  for(typename std::set<BBox<Type>*>::iterator it = unique.begin(); it != unique.end(); ++it){
    *vit = *it;
    vit++;
  }
}

template <class Type>
BBox<Type>::BBox()
{
  octree = NULL;
  level = 0;
  count = 0;
  nBelow = 0;

  Down = NULL;
  max[0] = max[1] = max[2] = 0.0;
  min[0] = min[1] = min[2] = 0.0;
 
  return;
}

template <class Type>
BBox<Type>::~BBox(){

  //delete the children of this object
  delete [] Down;

  return;
}

template <class Type>
void BBox<Type>::Split()
{
  Int i;
  nBelow = 8;
  Down = new BBox<Type>*[8];
 
  for(i = 0; i < 8; i++){
    Down[i] = new BBox<Type>;
    Down[i]->octree = octree;
    Down[i]->level = this->level + 1;
  }
  //adjust static depth
  if((this->level +1) > octree->deepestLevel) octree->deepestLevel = this->level+1;

  //compute the max and min extents of each level
  Type mid[3];
  for(i = 0; i < 3; i++){
    mid[i] = (this->max[i] - this->min[i])/2.0 + this->min[i];
  }

  //----------------------------------------------------------//
  //   y      OCTREE DIVISION PATTERNS
  //   |   z
  //   |  /                 FRONT:  2---3    BACK: 6---7
  //   | /                          |   |          |   |
  //   |/______x                    0---1          4---5
  //   
  //

  //Child 0
  //front plane - lower left corner
  Down[0]->min[0] = this->min[0];
  Down[0]->min[1] = this->min[1];
  Down[0]->min[2] = this->min[2];
  Down[0]->max[0] = mid[0];
  Down[0]->max[1] = mid[1];
  Down[0]->max[2] = mid[2];
  //Child 1
  //front plane - lower right corner
  Down[1]->min[0] = mid[0];
  Down[1]->min[1] = this->min[1];
  Down[1]->min[2] = this->min[2];
  Down[1]->max[0] = this->max[0];
  Down[1]->max[1] = mid[1];
  Down[1]->max[2] = mid[2];
  //Child 2
  //front plane - upper left corner
  Down[2]->min[0] = this->min[0];
  Down[2]->min[1] = mid[1];
  Down[2]->min[2] = this->min[2];
  Down[2]->max[0] = mid[0];
  Down[2]->max[1] = this->max[1];
  Down[2]->max[2] = mid[2];
  //Child 3
  //front plane - upper right corner
  Down[3]->min[0] = mid[0];
  Down[3]->min[1] = mid[1];
  Down[3]->min[2] = this->min[2];
  Down[3]->max[0] = this->max[0];
  Down[3]->max[1] = this->max[1];
  Down[3]->max[2] = mid[2];
  //Child 4
  //back plane - lower left corner
  Down[4]->min[0] = this->min[0];
  Down[4]->min[1] = this->min[1];
  Down[4]->min[2] = mid[2];
  Down[4]->max[0] = mid[0];
  Down[4]->max[1] = mid[1];
  Down[4]->max[2] = this->max[2];
  //Child 5
  //back plane - lower right corner
  Down[5]->min[0] = mid[0];
  Down[5]->min[1] = this->min[1];
  Down[5]->min[2] = mid[2];
  Down[5]->max[0] = this->max[0];
  Down[5]->max[1] = mid[1];
  Down[5]->max[2] = this->max[2];
  //Child 6
  //back plane - upper left corner
  Down[6]->min[0] = this->min[0];
  Down[6]->min[1] = mid[1];
  Down[6]->min[2] = mid[2];
  Down[6]->max[0] = mid[0];
  Down[6]->max[1] = this->max[1];
  Down[6]->max[2] = this->max[2];
  //Child 7
  //back plane - upper right corner
  Down[7]->min[0] = mid[0];
  Down[7]->min[1] = mid[1];
  Down[7]->min[2] = mid[2];
  Down[7]->max[0] = this->max[0];
  Down[7]->max[1] = this->max[1];
  Down[7]->max[2] = this->max[2];

  for(i = 0; i < 8; i++){
    Down[i]->BuildList(this);
  }

  //now we check to see if we've created any zero occupancy boxes,
  //if we have prune them
  std::vector<BBox<Type>*> tmp;
  for(i = 0; i < 8; i++){
    if(Down[i]->count == 0){
      delete Down[i];
      nBelow--;
    }
    else{
      tmp.push_back(Down[i]);
    }
  }
  
  delete [] Down;
  Down = new BBox<Type>*[8];
  for(Int i = 0; i < tmp.size(); i++){
    Down[i] = tmp[i];
  }

  //add the children to our flatlist
  for(i = 0; i < nBelow; i++){
    octree->flatList.push_back(Down[i]);
  }

  //free up the memory held by this bbox list, we don't need it after the split
  std::vector<Int> tempVect;
  list.swap(tempVect);

  return;
}

template <class Type>
Bool BBox<Type>::IsSplit()
{
  if(Down == NULL) return false;
  return true;
}

template <class Type>
void BBox<Type>::BuildList(BBox<Type>* Up)
{
  Int i;

  //guess at 1/8th the size of the parent level
  Int memSize = Up->count / 8 + 2;

  list.reserve(memSize);
  //only search over the nodes which are known to be in the parent list
  count = 0;
  Type* xyz = octree->xyz;
  for(i = 0; i < Up->count; i++){
    Int node = Up->list[i];
    if(real(xyz[node*3 + 0]) <= real(max[0]) && real(xyz[node*3 + 0]) >= real(min[0])){
      if(real(xyz[node*3 + 1]) <= real(max[1]) && real(xyz[node*3 + 1]) >= real(min[1])){
	if(real(xyz[node*3 + 2]) <= real(max[2]) && real(xyz[node*3 + 2]) >= real(min[2])){
	  //if the node is within the extents, add it to the list
	  list.push_back(node);
	  count++;
	}
      }
    }
  }
  
  return;
}

template <class Type>
BBox<Type>* BBox<Type>::GetLeaf(Type* xyz, Type* dist)
{
  Int i;
  Bool bottom;
  BBox<Type>* candidate = this;
  BBox<Type>* candidateNew = this;
  Type minDist, newDist;
  
  //should only be called from the root level
  if(level != 0){
    std::cerr << "WARNING: Cannot traverse tree starting from a child" << std::endl;
    return candidate;
  }
  
  //find minimum distance (non-empty) voxel and descend into it
  do{
    candidate = candidateNew;
    if(candidate->IsSplit()){
      minDist = 1.0e99;
      for(i = 0; i < candidate->nBelow; i++){
	newDist = candidate->Down[i]->GetVoxelDistance(xyz);
	if((real(newDist) < real(minDist)) && (candidate->Down[i]->count != 0)){
	  minDist = newDist;
	  candidateNew = candidate->Down[i];
	}
      }
      bottom = false;
    }
    else{
      bottom = true;
    }
  }while(!bottom);
  
  *dist = minDist;
  return candidate;
}

template <class Type>
Type BBox<Type>::GetNodeDistance(Type* xyz, Int& nodeId)
{
  Int i;
  Type distMin;
  Type minLeafDist;
  Type leafDist;
  Bool done;

  BBox<Type>* leaf = GetLeaf(xyz, &leafDist);
  BBox<Type>* candidate;

  //get the possible neighbors of this leaf
  std::vector<BBox<Type>*> neighbors;
  std::set<BBox<Type>*> attempted;
  nodeId = -1;

  //check for zero occupancy, just in case
  if(leaf->count == 0){
    std::cerr << "WARNING: a zero occupancy leaf was returned from GetLeaf()" << std::endl;
    nodeId = -1;
    return 1.0e99;
  }
  //loop over nodes in leaf and compute minimum distance
  distMin = 1.0e99;

  for(i = 0; i < leaf->list.size(); i++){
    Int node = leaf->list[i];
    Type dist = Distance(&octree->xyz[node*3 + 0], xyz);
    if(real(dist) < real(distMin)){
      distMin = dist;
      nodeId = node;
    }
  }
  leaf->GetNeighbors(neighbors, xyz, distMin);
  do{
    done = true;
    //now check each of the neighboring leafs to see if the point
    //of interest is closer to any of them than the nearest node
    minLeafDist = 1.0e99;
    for(i = neighbors.size()-1; i >= 0; i--){
      candidate = neighbors[i];
      neighbors.pop_back();
      //check for set insertion success, if we haven't tried it, give it a shot
      if(attempted.insert(candidate).second){
	if(real(candidate->GetVoxelDistance(xyz)) < real(distMin)){
	  leaf = candidate;
	  //signal to continue checking
	  done = false;
	  //break out of the loop, we have a new candidate to check now
	  break;
	}
      }
    }
    for(i = 0; i < leaf->list.size(); i++){
      Int node = leaf->list[i];
      Type dist = Distance(&octree->xyz[node*3 + 0], xyz);
      if(real(dist) < real(distMin)){
	distMin = dist;
	nodeId = node;
      }
    }
  }while(!done);

  return distMin;
}

template <class Type>
Type BBox<Type>::GetVoxelDistanceCentroid(const Type* xyz)
{
  Type x = 0.5*(min[0] + max[0]);
  Type y = 0.5*(min[1] + max[1]);
  Type z = 0.5*(min[2] + max[2]);

  Type dx = xyz[0] - x;
  Type dy = xyz[1] - y;
  Type dz = xyz[2] - z;

  Type dist = sqrt(dx*dx + dy*dy + dz*dz);

  return dist;
  
}

template <class Type>
Type BBox<Type>::GetVoxelDistance(const Type* xyz)
{
  Type dist;
  
  Type x = xyz[0];
  Type y = xyz[1];
  Type z = xyz[2];

  //find the edge of the box our point is closest to
  if(real(x) <= real(min[0])){
    x = min[0];
  }
  if(real(x) >= real(max[0])){
    x = max[0];
  }
  if(real(y) <= real(min[1])){
    y = min[1];
  }
  if(real(y) >= real(max[1])){
    y = max[1];
  }
  if(real(z) <= real(min[2])){
    z = min[2];
  }
  if(real(z) >= real(max[2])){
    z = max[2];
  }
  
  Type dx = xyz[0] - x;
  Type dy = xyz[1] - y;
  Type dz = xyz[2] - z;

  dist = sqrt(dx*dx + dy*dy + dz*dz);

  return dist;
}

template <class Type>
void BBox<Type>::GetExtents(Type* box)
{
  box[0] = min[0];
  box[1] = max[0];
  box[2] = min[1];
  box[3] = max[1];
  box[4] = min[2];
  box[5] = max[2];
}

template <class Type>
void BBox<Type>::Contract()
{
  //if this box is not at the bottom of the list, contract based on contained boxes
  if(Down != NULL){
    Type bounds[6];
    
    for(Int i = 0; i < 3; i++){
      max[i] = -1.0e99;
      min[i] = 1.0e99;
    }

    for(Int i = 0; i < nBelow; i++){
      BBox<Type>& box = *Down[i];
      //contract past boxes with zero occupancy, we want minimum volume fill
      if(box.count != 0){
	box.GetExtents(bounds);
	for(Int i = 0; i < 3; i++){
	  if(real(bounds[i*2 + 0]) < real(min[i])) min[i] = bounds[i*2 + 0];
	  if(real(bounds[i*2 + 1]) > real(max[i])) max[i] = bounds[i*2 + 1];
	}
      }
    }
    return;
  }
  if(count == 0){
    //if we contain no points contract to a singular value at the leaf centroid,
    //this allows upper level boxes to contract past this one
    for(Int i = 0; i < 3; i++){
      min[i] = max[i] = 0.5*(max[i]+min[i]);
    }
    return;
  }
  Type max[3];
  Type min[3];
  for(Int i = 0; i < 3; i++){
    max[i] = -1.0e99;
    min[i] = 1.0e99;
  }
  Type* xyz = octree->xyz;
  for(Int i = 0; i < count; i++){
    Int node = list[i];
    Type x = xyz[node*3 + 0];
    Type y = xyz[node*3 + 1];
    Type z = xyz[node*3 + 2];
    if(real(x) < real(min[0])) min[0] = x;
    if(real(x) > real(max[0])) max[0] = x;
    if(real(y) < real(min[1])) min[1] = y;
    if(real(y) > real(max[1])) max[1] = y;
    if(real(z) < real(min[2])) min[2] = z;
    if(real(z) > real(max[2])) max[2] = z;
  }
}
