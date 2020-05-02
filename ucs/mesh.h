#ifndef MESH_H__
#define MESH_H__

#include "general.h"
#include "uns_base.h"
#include "etypes.h"
#include "elementClass.h"
#include "statflags.h"
#include "endian_util.h"
#include "strings_util.h"
#include "tinyxml.h"
#include "parallel.h"
#include "h5layer.h"
#include "dataInfo.h"
#include "solutionField.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <mpi.h>
//used for Cuthill-McKee -- STL
#include <deque> 
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

#ifdef _HAS_CGNS
#include <cgnslib.h>
#endif

#include "geometry.h"
#include "mem_util.h"
#include "macros.h"
#include "oddsNends.h"
#include "exceptions.h"

void TranslateWinding(Int* nodes, Int translation[6][8], Int num_nodes, Int etype, Int to_other_format);

//forward declarations
template <class Type> class PObj;
template <class Type> class SolutionField;

template <class Type> 
class Mesh
{
  //this little bit of weirdness brought to you by c++
  //when c++ is your hammer everything looks like your thumb...
  //truly this just allows complex and double typed classes to use private members in copy constructors
  template <class Type3>
    friend class Mesh;
  
 public:

  //*************************************************************************************
  //     MEMBER FUNCTIONS
  //*************************************************************************************
  
  Mesh();
  ~Mesh();

  template <class Type2>
  Mesh(const Mesh<Type2>& meshToCopy);

  void SetParallelPointer(PObj<Type>* p_passed); //pass in parallel pointer for the mesh to use
  Int MemInitMesh();
  
  //This uses list ordering to reorder the list in c2n
  //By modifying c2n we can control how all other maps are built
  //to minimize sparse matrix bandwidth
  Int ReorderC2nMap();

  Int BuildMaps();   //MASTER CALLING FUNCTION - builds all conn. maps
                     //all other mapping functions order dependent - privatized
  
  Int BuildMapsDecomp(); //build maps only up to psp() primitives -- for DECOMP
  void CleanMapsDecomp(); //clean memory for maps up to psp() primitives -- for DECOMP

  void AllocateSolutionMemory();
  void KillMapMemory();       //deletes all memory associated with mapsInit()
  void KillSolutionMemory();  //deletes all memory associated with solution allocation
  void KillEdgesMemory();     //deletes all memory associated with edges in mesh

  Int CalcMetrics(); //MASTER CALLING FUNCTION - calcs all metrics, and checks
                     //all other metrics functions privatized 
  Int ReorderMeshCuthillMcKee(Int reverse); //Reorder mesh using Cuthill-McKee algorithm
  void ScaleMesh(Type scaleFactor);   

  void AppendNodes(Int numNewNodes, Type* xyz_new);
  
  //allocates array and returns ids of points with desired factag
  Int FindPointsWithFactag(Int** pts, Int factag);
  //allocates array and returns ids of elements with desired factag
  Int FindSurfaceElementsWithFactag(std::vector<Int>& elementIds, Int factag);
  //returns the average sizing of elements on a given factag
  Type ComputeElementSizingAverageOnFactag(Int factag);

  Bool IsScaled() const {return scaled;};
  Bool IsReordered() const {return reordered;};
  Bool IsElemParallel(Int gelem);
  Int IsInteriorNode(Int n) const;   //interior to mesh
  Int IsGhostNode(Int n) const;      //parallel update nodes
  Int IsBoundaryNode(Int n) const;   //boundary condition nodes
  Int GetNumElemNodes(Int type) const {return mnode[type];};
  Int GetNumElem() const;
  Int GetNumElem(Int type) const {return nelem[type];};
  Int GetNumParallelNodes() const {return gnode;};
  Int GetNumGlobalNodes() const {return gnnode;};
  Int GetNumNodes() const {return nnode;};
  Int GetNumBoundaryNodes() const {return nbnode;};
  Int GetNumGlobalElem() const {return gnelem;};
  Int GetNumLocalElem() const {return lnelem;};
  Int GetNumEdges() const {return nedge;};
  Int GetNumParallelEdges() const {return ngedge;};
  Int GetNumBoundaryEdges() const {return nbedge;};
  Type GetVolumeTotal() const;
  Type const * GetVolume() const;       //returns vector of cv volumes
  Type const * GetVolumeOld() const;    //returns vector of nm1 cv volumes
  Type const * GetVolumeOldM1() const;  //returns vector of nm2 cv volumes
  Type const * GetNodeCoords() const;
  Type const * GetNodeCoordsBase() const;
  Int GetMaximumFactag();
  void GetCoordsForPoints(Type* rxyz, Int* pts, Int n); 
  Type GetMinX() const {return extentsMin[0];};
  Type GetMaxX() const {return extentsMax[1];};
  Type GetMinY() const {return extentsMin[2];};
  Type GetMaxY() const {return extentsMax[0];};
  Type GetMinZ() const {return extentsMin[1];};
  Type GetMaxZ() const {return extentsMax[2];};
  //give local node number, returns owning process and the localid of that node 
  //on the owning process... owning Process and localId return negative if in error
  void GetParallelTuple(Int node, Int& owningProcess, Int& localId) const;
  //computes the normal vector to a node based on neighbor faces
  void GetNodeNeighborhoodNormal(Int ptid, std::vector<Int> excludedFactags, std::vector<Type>& normal) const; 
  
  void SetNumNodes(Int nnode){this->nnode = nnode;};
  void SetMeshScaled(){scaled  = true;};
  void SetMeshReordered(){reordered = true;};
  Type* GetNodeCoords();
  Int* ElspBegin(Int ptid); //returns beginning access to elements surrounding a point
  Int* ElspEnd(Int ptid);   //returns ending access to elements surrounding a point
  Int* SelspBegin(Int ptid); //returns beginning access to elements surrounding a surface point
  Int* SelspEnd(Int ptid);   //returns ending access to elements surrounding a surface point
  std::string GetStringFromElemType(Int etype); //returns a string stating the type of element from numerical type
  const Int* ElspBegin(Int ptid) const; //returns beginning access to elements surrounding a point
  const Int* ElspEnd(Int ptid) const;   //returns ending access to elements surrounding a point
  const Int* SelspBegin(Int ptid) const; //returns beginning access to elements surrounding a surface point
  const Int* SelspEnd(Int ptid) const;   //returns ending access to elements surrounding a surface point
  void GetInteriorFaces(std::vector<Element<Type>*>& newlist, std::vector<IntInt>& orientation);
  void GetBoundaryFaces(std::vector<Element<Type>*>& newlist, std::vector<Int>& volumeNeighbor,
			Int targetFactag);
  
  Int ReadPartedMesh(std::string casename);
  Int ReadUGRID_Ascii(std::string filename);
  Int ReadCRUNCH_Ascii(std::string filename);
  Int ReadSU2_Ascii(std::string filename);
  Int ReadGMSH_Master(std::string filename);
  Int ReadGMSH2_Ascii(std::string filename);
  Int ReadGMSH40_Ascii(std::string filename);
  Int ReadGMSH41_Ascii(std::string filename);
  void EliminateOrphanedNodes();
  void UpdateElementCounts(); // updates the element type counters
#ifdef _HAS_CGNS
  Int ReadCGNS(std::string filename); // master function calling utilities below
  Int GetCGNSSizes(std::string filename, cgsize_t** isize);
  void CGNSreadCoordElem(std::string filename, cgsize_t**isize);
  int CGNSgetBCNum(char* name);
  void CGNSreadBCConditions(char *name, int **bc);
  int CGNSgetBCIndNum(char *name, int ib);
#endif
  Int WriteSTL_Ascii(std::string casename, Int selectFactag);
  Int WriteFluentCase_Ascii(std::string casename);
  Int WriteCRUNCH_Ascii(std::string casename);
  Int WriteVTK_Ascii(std::string casename, std::vector<SolutionField<Type>*>& fields);
  Int WriteVTK_Binary(std::string casename, std::vector<SolutionField<Type>*>& fields);
  Int WriteXMLVTK_Binary(std::string casename, std::vector<SolutionField<Type>*>& fields);
  Int WriteGridXDMF(PObj<Type> &p, std::string filebase, std::string meshname);
  Int WriteUGRID_Ascii(std::string casename);
  //this call will write a decomp file as appropriate for 
  //this portion of the mesh in its current state
  Int WriteParallelMesh(std::string casename);
  //this call writes current coords in a timestep-# directory for mesh movement
  Int WriteCurrentCoords(std::string casename, Int timestep);

  std::string i2h(Int a){
    std::stringstream ss;
    ss << std::hex << a;
    return ss.str();
  };
  
  //*************************************************************************************
  //     MEMBER VARIABLES
  //*************************************************************************************
  
  
  PObj<Type>* p;  //parallel comm object... we need one of thes for each mesh we use

  std::vector<Element<Type>*> elementList; //list of mesh elements
  Edges<Type>* edges;                      //list of Edge objects
  HalfEdges<Type>* bedges;                 //list of HalfEdge objects (ghost edges then boundary edges)
  
  //parallel utility maps
  Int* gNodeOwner;   //(gnode long) lists proc id which own ghost nodes
  Int* gNodeLocalId; //(gnode long) lists location of node in local owner list
                     //this is the node's location on the process which owns it
                     //this list is ordered by process... ie. all of process 1's 
                     //    nodes are listed together, etc.

  //solution vectors
  Type* cg;    //center of gravity
  Type* nv;    //node velocity

  //edge coefficients for LSQ gradients
  Type* s;
  Type* sw;

  //solution ordering
  Int* ordering;
  
  Int nbface;    //number of surface elements (TRI/QUAD)
  Int nfactags;  //number of unique face tags in list factag

  Type* xyz;      //coordinates of each vertex time current
  Type* xyz_base; //coordinates of each vertex at t=0, needed for movement routines
  Type* xyzold;   //position time n 
  Type* xyzoldm1; //position time n-1

  Type* vol;      //volumes of node based CVs
  Type* volold;   //volumes at n   (for GCL)
  Type* vololdm1; //volumes at n-1 (for GCL)
  
  Int* psp;      //point surrounding point map
  Int* ipsp;     //index into CRS of psp
  Int* besp;     //boundary halfedges surrounding point
  Int* ibesp;    //index into CRS of besp

  Int nvnodes;   //number of nodes not lying on a surface
  Int* vnodes;   //simple list of nodes not on a surface (in field only)
  Int nsnodes;   //number of nodes lying on a surface
  Int* snodes;   //simple list of nodes directly on a surface

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

private:

  //*************************************************************************************
  //     PRIVATE MEMBER FUNCTIONS
  //*************************************************************************************

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
  Int GetSharedFaceNodes(Element<Type>& elem1, Element<Type>& elem2, Int* knodes, Int& order);
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


  //*************************************************************************************
  //     PRIVATE MEMBER VARIABLES
  //*************************************************************************************
  
  //number of nodes in each element type
  Int* mnode;
  Type* extentsMax; //maximum xyz coords
  Type* extentsMin; //minimum xyz coords
  
  //list of booleans which indicates what memory has been alloc'd
  Bool meshInit;
  Bool mapsInit;
  Bool mapsInitPsp;
  Bool edgeInit;
  Bool metricsCalc;
  Bool solMemAlloc;
  //flag that is set on mesh read which states a mesh is previously
  //scaled.. we wont' do it again
  Bool scaled;
  //flag that is set on mesh read which states a mesh is previously
  //reordered... we won't do it again
  Bool reordered;

 protected:

  //*************************************************************************************
  //     PROTECTED MEMBER FUNCTIONS
  //*************************************************************************************


  //*************************************************************************************
  //     PROTECTED MEMBER VARIABLES
  //*************************************************************************************
  
  Int* esp;      //edge surrounding point
  Int* iesp;     //index into CRS of esp

  Int* el2el;   //volume element to element connectivity (includes surf. elements)
  Int* iel2el;  //index into CRS of vel2el

  Int* sel2el;   //surface element to element connectivity (excludes vol. elements)
  Int* isel2el;  //index into CRS of sel2el

  Int* elsp;     //elements surrounding point map (includes surf. elemente)
  Int* ielsp;    //index into CRS of elsp
  Int* selsp;    //elements surrounding point map (exludes vol. elements)
  Int* iselsp;   //index into CRS of selsp
};

//include implementations
#include "mesh.tcc"

#endif


