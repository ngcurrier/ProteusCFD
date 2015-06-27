#ifndef MESH_H__
#define MESH_H__

#include "general.h"
#include "io.h"
#include "uns_base.h"
#include "elementClass.h"
#include "statflags.h"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>

//forward declarations
template <class Type> class PObj;

template <class Type> 
class Mesh
{
  
public:

  Mesh();
  ~Mesh();

  template <class Type2>
  Mesh(const Mesh<Type2>& meshToCopy);

  PObj<Type>* p;  //parallel comm object... we need one of thes for each mesh we use
  void SetParallelPointer(PObj<Type>* p_passed); //pass in parallel pointer for the mesh to use
  Int ReadPartedMesh(std::string casename);

  Int MemInitMesh();
  
  Int ReorderC2nMap();//This uses list ordering to reorder the list in c2n
                     //By modifying c2n we can control how all other maps are built
                     //to minimize sparse matrix bandwidth

  Int BuildMaps();   //MASTER CALLING FUNCTION - builds all conn. maps
                     //all other mapping functions order dependent - privatized
  
  Int BuildMapsDecomp(); //build maps only up to psp() primitives -- for DECOMP
  void CleanMapsDecomp(); //clean memory for maps up to psp() primitives -- for DECOMP

  void KillMapMemory(); //deletes all memory associated with mapsInit()
  void KillSolutionMemory(); //deletes all memory associated with solution allocation
  void KillEdgesMemory();   //deletes all memory associated with edges in mesh

  //adaptation routines
  void RefineNodeBased(Bool* refineNodeList);
  void RebuildRefinedMesh(Type* xyzNew, Int xyzNewCount, Element<Type>** newElements, Int newElemCount, Int* deadElemList, Int deadElemCount);

  Bool IsElemParallel(Int gelem);

  //give local node number, returns owning process and the localid of that node 
  //on the owning process... owning Process and localId return negative if in error
  void GetParallelTuple(Int node, Int* owningProcess, Int* localId);

  Int IsInteriorNode(Int n) const;   //interior to mesh
  Int IsGhostNode(Int n) const;      //parallel update nodes
  Int IsBoundaryNode(Int n) const;   //boundary condition nodes

  Int CalcMetrics(); //MASTER CALLING FUNCTION - calcs all metrics, and checks
                     //all other metrics functions privatized 

  Int ReorderMeshCuthillMcKee(Int reverse); //Reorder mesh using Cuthill-McKee algorithm
  
  Int AllocateSolutionMemory();

  Int NodeOrigToLocal(Int origId);

  //allocates array and returns ids of points with desired factag
  Int FindPointsWithFactag(Int** pts, Int factag);
  void GetCoordsForPoints(Type* rxyz, Int* pts, Int n); 

  //list of elements in the mesh
  std::vector<Element<Type>*> elementList;

  //this call will write a decomp file as appropriate for 
  //this portion of the mesh in its current state
  Int WriteParallelMesh(std::string casename);
  //this call writes current coords in a timestep-# directory for mesh movement
  Int WriteCurrentCoords(std::string casename, Int timestep);

  //scales mesh by a value
  void ScaleMesh(Type scaleFactor);

  //flag that is set on mesh read which states a mesh is previously
  //reordered... we won't do it again
  Bool reordered;

  //flag that is set on mesh read which states a mesh is previously
  //scaled.. we wont' do it again
  Bool scaled;

  //number of nodes in each element type
  Int* mnode;
  
  Int gnnode;       //number of global nodes
  Int nnode;        //number of nodes in mesh
  Int gnode;        //number of ghost nodes in mesh (parallel only)
  Int nbnode;       //nodes numbered >= (nnode+gnode) & < (nnode+gnode+nbnode) are phantom/ghost nodes for bedges
  Int* nelem;       //number of each type of element in mesh -- nelem[Tri], nelem[Quad], ...
  Int gnelem;       //number of global elements 
  Int lnelem;       //number of elements in local mesh
  Int nedge;        //number of edges in mesh
  Int ngedge;       //number of half edges in mesh (ghost -- for parallel)
  Int nbedge;       //number of half edges in mesh (boundary)

  //parallel utility maps
  Int* gNodeOwner;   //(gnode long) lists proc id which own ghost nodes
  Int* gNodeLocalId; //(gnode long) lists location of node in local owner list
                     //this is the node's location on the process which owns it
                     //this list is ordered by process... ie. all of process 1's 
                     //    nodes are listed together, etc.

  //solution vectors
  Type* cg;    //center of gravity
  Type* nv;    //node velocity
  Type* xyzold;   //position time n 
  Type* xyzoldm1; //position time n-1

  //edge coefficients for LSQ gradients
  Type* s;
  Type* sw;

  //solution ordering
  Int* ordering;

  //number of surface elements (TRI/QUAD)
  Int bface;

  Int nfactags;  //number of unique face tags in list factag

  Type* extentsMax; //maximum xyz coords
  Type* extentsMin; //minimum xyz coords

  Type* xyz;      //coordinates of each vertex
  Type* xyz_base; //coordinates of each vertex at t=0, needed for movement routines
  Type* vol;      //volumes of node based CVs
  Type* volold;   //volumes at n   (for GCL)
  Type* vololdm1; //volumes at n-1 (for GCL)

  CVStat* cvstat; //control volume boolean status flags

  Int* elsp;     //elements surrounding point map (includes surf. elemente)
  Int* ielsp;    //index into CRS of elsp
  Int* selsp;    //elements surrounding point map (exludes vol. elements)
  Int* iselsp;   //index into CRS of selsp

  Int* psp;      //point surrounding point map
  Int* ipsp;     //index into CRS of psp

  Int* el2el;   //volume element to element connectivity (includes surf. elements)
  Int* iel2el;  //index into CRS of vel2el

  Int* sel2el;   //surface element to element connectivity (excludes vol. elements)
  Int* isel2el;  //index into CRS of sel2el

  Int nvnodes;   //number of nodes not lying on a surface
  Int* vnodes;   //simple list of nodes not on a surface (in field only)
  Int nsnodes;   //number of nodes lying on a surface
  Int* snodes;   //simple list of nodes directly on a surface

  Edges<Type>* edges;       //list of Edge objects
  HalfEdges<Type>* bedges;  //list of HalfEdge objects (ghost edges then boundary edges)
  
  Int* esp;      //edge surrounding point
  Int* iesp;     //index into CRS of esp

  Int* besp;     //boundary halfedges surrounding point
  Int* ibesp;    //index into CRS of besp

private:
  Bool meshInit;
  Bool mapsInit;
  Bool mapsInitPsp;
  Bool edgeInit;
  Bool metricsCalc;
  Bool solMemAlloc;
  
  //Called via public BuildMaps() function
  Int BuildElsp();  //element to surrounding point map
  Int BuildPsp();   //point to surrounding point map
  Int BuildEl2el(); //element to element map
  Int FixBoundaryWindings(); //rewind boundary elements so the normal point inward consistently
  Int BuildSel2el(); //surface element to element map (called only if b. windings sane)
  Int BuildEdges(); //edge objects which contain jacobians, etc.
  Int BuildEsp();    //edge to point map with parallel features

  //Internal Utility Functions
  Int EliminateRepeats(Int* array, Int n);
  Int CheckSharedFace(Element<Type>& elem1, Element<Type>& elem2, Int* knodes);
  Int CheckFaceLocation(Int** fli, Int* fl, Int etype, Int facenum, Int* knodes, Int nodes);
  Int CheckSharedEdge(Element<Type>& elem1, Element<Type>& elem2, Int* knodes);
  Int CheckEdgeLocationSurf(Int** eli, Int* el, Int etype, Int ednum, Int* knodes, Int nodes);
  Int EdgeInElem(Int n1, Int n2, Int jelem);
  Int CheckEdgeOrientation(Int** ieorder, Int* eorder, Int etype, Int numedges, 
			   Int* tested, Int& edgeid);
    
  //metric calculation routines
  Int CalcExtents();         //calculates extents of mesh... 
  Int CalcAreasVolumes();    //calculate area vectors for the edges, and volumes of CVs
  Int CheckForClosure();     //check for CV closure

};

//include implementations
#include "mesh.tcc"
#include "meshAdapt.tcc"

#endif


