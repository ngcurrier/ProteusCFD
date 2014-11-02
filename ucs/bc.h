#ifndef BC_H__
#define BC_H__

#include "general.h"
#include "kernel.h"
#include "bcobj.h"
#include "powerLaw.h"
#include "elementClass.h"
#include <string>
#include <vector>

//forward declarations
template <class Type> class Mesh;
template <class Type> class Param;
template <class Type> class EqnSet;

//class which stores data relevant
//to all present boundary conditions
template <class Type>
class BoundaryConditions
{
 public:
  BoundaryConditions();
  ~BoundaryConditions();

  Int ReadFile(std::string casename);
  void PrintBCs();

  //returns the numerical bc type id with input of the factag id
  Int GetBCId(Int factag);
  //returns the correct BCObj pointer with input of the factag id
  BCObj<Type>* GetBCObj(Int factag);

  Int InitInternals(EqnSet<Type>* eqnset); //this sets individual BCObj pointers internally, etc.

  //returns true if the node is connected to a face of type given
  template <class Type2>
  Bool IsNodeOnBC(Int nodeid, Mesh<Type2>* m, Int bcId);

  std::string* surface_names;   //plain text identifiers for surfaces
  Int* bc_applied;              //list which holds bc's applied to surfaces
  Int* bc_map;                  //maps bc number in file to lower order number in code
                                //this is in case number is non-consecutive i.e. 1, 3, 199, etc.

  BCObj<Type>* bcs;                   //list of reference states

  Int num_bcs;                  //number of bcs present in domain
  Int largest_bc_id;            //largest bc id present in bc file
  Int num_bodies;               //number of composite bodies
  Int largest_body_id;          //largest composite body id found

 private:
 
  void CountNumberBCsInFile(std::ifstream& fin); //sets number of BCs in file, return largest id found
  void CountBodiesInFile(std::ifstream& fin);  //sets number of composite bodies in file
  Int ParseLine(std::string& line);
  Int ParseBody(std::string& line);
  Int SetVarsFromLine(Int ObjId, std::string& line);  //parse variables from line into BCObj
  void AllocateBCObjects();     //allocate required bc objects - one per bc

};


//this calls out any jacobian modifications that need to be made due to BCs
template <class Type>
void Bkernel_BC_Jac_Modify(B_KERNEL_ARGS);

//this calls out any residual modifications that need to be made due to BCs
template <class Type>
void Bkernel_BC_Res_Modify(B_KERNEL_ARGS);

//this templating is getting a bit silly but we have to have it to support
//complex jacobians...
template <class Type, class Type2, class Type3>
void CalculateBoundaryVariables(EqnSet<Type>* eqnset, Mesh<Type3>* m, SolutionSpace<Type3>* space, 
				Type* QL, Type* QR, Type* Qinf, Type* avec, Int bcId, 
				Int eid, BCObj<Type2>* bcobj, Type3 vdotn, Type3* velw);

//this is the master call to update the boundary conditions
template <class Type>
void UpdateBCs(SolutionSpace<Type>* space);

//this is the kernel that actually updates the BCs
template <class Type>
void BC_Kernel(B_KERNEL_ARGS);

//function returns number of points strictly on symmetry planes... used for smoothing meshes
//which have been extruded in 2D.. allocates memory.. pass NULL pointer
template <class Type>
Int GetNodesOnSymmetryPlanes(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int** nodes);

//function returns number of points on a particular type of bc... used for smoothing meshes
//allocates memory.. pass NULL pointer
template <class Type>
Int GetNodesOnBCType(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int** nodes, 
		     Int BCType);

//function returns number of points which are on a bc with any design tag
//allocates memory.. pass NULL pointer
template <class Type>
Int GetNodesOnMovingBC(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int** nodes);

//function returns number of points on a particular design BC (beta)... used for computational
//design work.. allocates memory.. pass NULL pointer
template <class Type>
Int GetNodesOnMovingBC(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int beta, 
		       Int** nodes);

//function returns number of points on a particular bcid
//allocates memory pass NULL pointer
template <class Type>
Int GetNodesOnBCId(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int bcid, Int** nodes);

//function takes a node list and returns a unique list of surface elements which
//those nodes are connected to
template <class Type>
void GetSElemsFromNodes(const Mesh<Type>* m, BoundaryConditions<Real>* bc, Int bcid, 
			std::vector<Element<Type>*>& elementList);

//sets status flags for each cv based on boundary conditions
template <class Type>
void SetBCStatFlags(Mesh<Type>* m, BoundaryConditions<Real>* bc);

//include implementations
#include "bc.tcc"

#endif
