#include "geometry.h"
#include "elementClass.h"
#include "mem_util.h"
#include "etypes.h"
#include "macros.h"
#include "parallel.h"
#include "oddsNends.h"
#include "h5layer.h"
#include "exceptions.h"
#include <mpi.h>
//used for Cuthill-McKee -- STL
#include <deque> 
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

//used to check if a list is negative or not
//class must allow random access via []
template<class theType>
Int NegativeQueue(theType* A, Int n){
  Int i;
  for(i = 0; i < n; i++){
    if(A[i].a >= 0){
      return(0);
    }
  }
  return(1);
}

//default constructor
template <class Type> 
Mesh<Type>::Mesh()
{
  //set some default values
  nnode = 0;
  gnode = 0;
  nbnode = 0;
  nbedge = 0;
  ngedge = 0;

  //allocate space for element counts
  nelem = new Int[MAX_E_TYPES];
  mnode = new Int[MAX_E_TYPES];

  //allocate space for extents
  extentsMax = new Type[3];
  extentsMin = new Type[3];

  //set number of vertices to expect for each elem type
  mnode[TRI] = 3;
  mnode[QUAD] = 4;
  mnode[TET] = 4;
  mnode[PYRAMID] = 5;
  mnode[PRISM] = 6;
  mnode[HEX] = 8;

  meshInit = 0;
  mapsInit = 0;
  mapsInitPsp = 0;
  edgeInit = 0;
  metricsCalc = 0;
  
  solMemAlloc = 0;
  
  reordered = false;
  scaled = false;

  edges = NULL;
  bedges = NULL;
  isel2el = NULL;
  sel2el = NULL;
}

template <class Type>
void Mesh<Type>::SetParallelPointer(PObj<Type>* p_passed)
{
  p = p_passed;
  return;
}

template <class Type>
Int Mesh<Type>::MemInitMesh()
{
  if(nnode == 0){
    std::cerr << "Memory allocation undefined: Mesh::memInit()" << std::endl;
    return (-1);
  }

  xyz = new Type[3*nnode + 3*gnode];
  xyz_base = new Type[3*nnode + 3*gnode];
  vol = new Type[nnode];
  volold = new Type[nnode];
  vololdm1 = new Type[nnode];
  
  //initliaze cv status flag objects
  cvstat = new CVStat[nnode];

  //allocate parallel memory
  gNodeOwner = new Int[gnode];
  gNodeLocalId = new Int[gnode];

  //this must be set here incase mesh reordering is not use
  //i.e. this must be a valid list as it is used in the buildMaps() function
  ordering = new Int[nnode];
  for(Int i = 0; i < nnode; i++){
    ordering[i] = i;
  }

  meshInit = 1;
  return (0);
}

template <class Type>
Int Mesh<Type>::ReorderC2nMap()
{
  Int err = 0;
  Int i;
  Int nodeorig, nodenew;
  Type* xyztemp = new Type[nnode*3];
  Int* origNodes = new Int[nnode];

  //Pick old node numbering, remap it, and reassign to c2n mapping
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes;
    Int nnodes = element.GetNodes(&nodes);
    for(i = 0; i < nnodes; i++){
      nodeorig = nodes[i];
      //make sure we are only moving interior nodes
      if(IsInteriorNode(nodeorig)){
	nodenew = ordering[nodeorig];
	nodes[i] = nodenew;
      }
    }
  }

  //switch out the coordinates for the new node
  memcpy(xyztemp, xyz, sizeof(Type)*nnode*3);
  for(i = 0; i < nnode; i++){
    nodenew = ordering[i];
    memcpy(&xyztemp[nodenew*3], &xyz[3*i], sizeof(Type)*3);
  }
  memcpy(xyz, xyztemp, sizeof(Type)*nnode*3);
  memcpy(xyz_base, xyztemp, sizeof(Type)*nnode*3);

  reordered = true;

  delete [] xyztemp;
  delete [] origNodes;

  return err;
}

template <class Type>
Int Mesh<Type>::BuildMaps()
{
  int err = 0;
  err = BuildElsp();
  if(err != 0){
    return (err);
  }
  err = BuildPsp();
  if(err != 0){
    return (err);
  }
  err = BuildEl2el();
  if(err != 0){
    return (err);
  }
  err = FixBoundaryWindings();
  if(err != 0){
    return (err);
  }
  err = BuildSel2el();
  if(err != 0){
    return (err);
  }
  err = BuildEdges();
  if(err != 0){
    return (err);
  }
  err = BuildEsp();
  if(err != 0){
    return (err);
  }

  AllocateSolutionMemory();
  
  mapsInit = 1;
  mapsInitPsp = 1;

  return(err);
}

//builds maps only up to psp() necessary for METIS interface
template <class Type>
Int Mesh<Type>::BuildMapsDecomp()
{
  int err = 0;
  err = BuildElsp();
  if(err != 0){
    return (err);
  }
  err = BuildPsp();
  if(err != 0){
    return (err);
  }
  
  mapsInitPsp = 1;
  return(err);
}
//clears memory for maps built with BuildMapsDecomp()
template<class Type>
void Mesh<Type>::CleanMapsDecomp()
{
  if(mapsInitPsp == 1){
    delete [] selsp;
    delete [] iselsp;
    delete [] elsp;
    delete [] ielsp;
    delete [] psp;
    delete [] ipsp;
    delete [] vnodes;
    delete [] snodes;
  }
  mapsInitPsp = 0;
}

template <class Type>
void Mesh<Type>::KillMapMemory()
{
  if(mapsInit == 1){
    //TODO: more robust memory freeing here - could still die
    delete [] el2el;
    delete [] iel2el;
    delete [] sel2el;
    delete [] isel2el;
    delete [] esp;
    delete [] iesp;
    delete [] besp;
    delete [] ibesp;
  }
  CleanMapsDecomp();
  mapsInit = mapsInitPsp = false;
  return;
}

template <class Type>
void Mesh<Type>::KillSolutionMemory()
{
  if(solMemAlloc == 1){
    delete [] s;
    delete [] sw;
    delete [] nv;
    delete [] xyzold;
    delete [] xyzoldm1;
  }
  solMemAlloc = false;
  return;
}


template <class Type>
void Mesh<Type>::KillEdgesMemory()
{
  if(edgeInit == 1){
    delete [] edges;
    delete [] bedges;
    delete [] cg;
  }
  edgeInit = false;

}

template <class Type>
Int Mesh<Type>::BuildElsp()
{
  Int node, indx, indxs, gelem;

  ielsp = NULL;
  elsp = NULL;
  iselsp = NULL;
  selsp = NULL;

  std::cout << "MESH UTILITY: Building element to surrounding point map... " << std::endl;

  ielsp = new Int[nnode+gnode+1];
  iselsp = new Int[nnode+gnode+1];

  //init CRS pointer
  for(Int i = 0; i < nnode+gnode+1; i++){
    ielsp[i] = 0;
    iselsp[i] = 0;
  }

  //fill in CRS pointer
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes;
    Int nnodes = element.GetNodes(&nodes);
    Int type = element.GetType();
    for(Int j = 0; j < nnodes; ++j){
      node = nodes[j];
      ielsp[node+1]++;
      if(type == TRI || type == QUAD){
	iselsp[node+1]++;
      }
    }
  }

  //setup offsets
  for(Int i = 1; i <= nnode+gnode; i++){
    ielsp[i] += ielsp[i-1];
    iselsp[i] += iselsp[i-1];
  } 

  //allocate storage for elsp map
  elsp = new Int[ielsp[nnode+gnode]];
  selsp = new Int[iselsp[nnode+gnode]];


  if(ielsp == NULL || elsp == NULL || iselsp == NULL || selsp == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildElsp()" << std::endl;
    return (-1);
  }
  
  //fill in map
  gelem = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes;
    Int nnodes = element.GetNodes(&nodes);
    Int type = element.GetType();
    for(Int j = 0; j < nnodes; ++j){
      node = nodes[j];
      indx = ielsp[node];
      elsp[indx] = gelem;
      ielsp[node]++;
      if(type == TRI || type == QUAD){
	indxs = iselsp[node];
	selsp[indxs] = gelem;
	iselsp[node]++;
      }
    }
    gelem++;
  }
  
  //shift indices in CRS back
  for(Int i = nnode+gnode; i > 0; i--){
    ielsp[i] = ielsp[i-1];
    iselsp[i] = iselsp[i-1];
  }
  ielsp[0] = 0;
  iselsp[0] = 0;

  //build list of nodes not lying on a surface
  //just loop through our surface elem surrounding point map and count
  nvnodes = 0;
  nsnodes = 0;
  for(Int i = 0; i < nnode+gnode; i++){
    if((iselsp[i+1] - iselsp[i]) == 0){
      nvnodes++;
    }
    else{
      nsnodes++;
    }
  }
  vnodes = new Int[nvnodes];
  snodes = new Int[nsnodes];
  Int j = 0;
  Int k = 0;
  for(Int i = 0; i < nnode+gnode; i++){
    if((iselsp[i+1] - iselsp[i]) == 0){
      vnodes[j] = i; 
      j++;
    }
    else{
      snodes[k] = i;
      k++;
    }
  }


  return (0);
}

template <class Type>
Int Mesh<Type>::BuildPsp()
{
  Int i,j,k;
  Int count1, count2;
  Int memsize, memjump, memcount, memmult;
  Int indx1, indx2;
  Int indx3, indx4;
  //element identifiers
  Int gelem, jnode;
  //local node number in element
  Int lnode = -1;

  //list of flags to perform sanity checks on connecting list
  Bool* checked;

  //set default allocation for psp and amount to 
  //add when resize takes place -- TUNABLE PARAMETERS
  memsize = 8*(nnode+gnode);
  memjump = (Int)real((1.5*(Type)(nnode+gnode)));
  memmult = 2;
  memcount = 1;

  //CRS storage for element point neighbor lookups 86 values
  //nodes connected to a particular point in a node
  Int cl[] = 
    {1, 2, 2, 0, 0, 1, //Tri 
     1, 3, 2, 0, 3, 1, 0, 2,  //Quad 
     1, 2, 3, 2, 0, 3, 0, 1, 3, 0, 1, 2,  //Tet 
     1, 3, 4, 2, 0, 4, 3, 1, 4, 0, 2, 4, 0, 1, 2, 3, //Pyramid 
     1, 2, 3, 2, 0, 4, 0, 1, 5, 4, 5, 0, 5, 3, 1, 3, 4, 2,     //Prism
     1, 3, 4, 2, 0, 5, 3, 1, 6, 0, 2, 7, 5, 7, 0, 6, 4, 1, 7, 5, 2, 4, 6, 3  //Hex
    };

  Int** cli = NULL;
  Int* cli_d = NULL;

  ipsp = NULL;
  psp = NULL;

  std::cout << "MESH UTILITY: Building point to surrounding point map... " << std::endl;

  cli = new Int*[6];
  cli_d = new Int[31];

  ipsp = new Int[nnode+gnode+1];
  ipsp[0] = 0;
  psp = new Int[memsize];

  if(cli == NULL || cli_d == NULL || ipsp == NULL || psp == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildPsp()" << std::endl;
    return (-1);
  }
  
  //setup lookup table by hand
  i = 0;
  cli[TRI] = &cli_d[i];
  i += mnode[TRI];
  cli[QUAD] = &cli_d[i];
  i += mnode[QUAD];
  cli[TET] = &cli_d[i];
  i += mnode[TET];
  cli[PYRAMID] = &cli_d[i];
  i += mnode[PYRAMID];
  cli[PRISM] = &cli_d[i];
  i += mnode[PRISM];
  cli[HEX] = &cli_d[i];

  //each node in an element requires two pointers...
  //one for the beginning of the connecting element
  //lookup and another for the end.. due to pyramids
  cli[TRI][0] = 0; 
  cli[TRI][1] = 2;
  cli[TRI][2] = 4;
  cli[QUAD][0] = 6;
  cli[QUAD][1] = 8;
  cli[QUAD][2] = 10;
  cli[QUAD][3] = 12;
  cli[TET][0] = 14;
  cli[TET][1] = 17;
  cli[TET][2] = 20;
  cli[TET][3] = 23;
  cli[PYRAMID][0] = 26;
  cli[PYRAMID][1] = 29;
  cli[PYRAMID][2] = 32;
  cli[PYRAMID][3] = 35;
  cli[PYRAMID][4] = 38;
  cli[PRISM][0] = 42;
  cli[PRISM][1] = 45;
  cli[PRISM][2] = 48;
  cli[PRISM][3] = 51;
  cli[PRISM][4] = 54;
  cli[PRISM][5] = 57;
  cli[HEX][0] = 60;
  cli[HEX][1] = 63;
  cli[HEX][2] = 66;
  cli[HEX][3] = 69;
  cli[HEX][4] = 72;
  cli[HEX][5] = 75;
  cli[HEX][6] = 78;
  cli[HEX][7] = 81;
  //Hex only has 8 nodes but we have to set the endpoint value...
  cli[HEX][8] = 84;

  count2 = 0;
  for(i = 0; i < nnode+gnode; i++){
    count1 = ipsp[i];
    indx1 = ielsp[i];
    indx2 = ielsp[i+1];
    for(j = indx1; j < indx2; j++){
      gelem = elsp[j];
      Element<Type>& element = *elementList[gelem];
      Int e = element.GetType();
      Int* nodes; 
      Int nnodes = element.GetNodes(&nodes);
      //find location of node i in ielem
      for(k = 0; k < nnodes; k++){
	if(i == nodes[k]){
	  lnode = k;
	  break;
	}
      }
      if(k == nnodes){
	std::cerr << "WARNING: In BuildPsp[] node " << i << " not found in element " << gelem << std::endl;
      }
      //find nodes connected to current node in own element
      //using connection lookup table cl, cli
      if((cli[e][lnode+1] - cli[e][lnode])+count2 >= memsize){
	MemResize(&psp, memsize, memsize+(memmult*memcount*memjump));
	memsize += (memmult*memcount*memjump);
	memcount++;
      }
      for(k = cli[e][lnode]; k < cli[e][lnode+1]; k++){
	jnode = nodes[cl[k]];
	psp[count2] = jnode;
	count2++;
      }
    }
    //eliminate duplicates from the list
    count2 = EliminateRepeats(&psp[count1], count2-count1);
    count2 += count1;
    ipsp[i+1] = count2;
  }

  //resize psp to be just big enough
  //using this function is a mem copy
  //and is likely more costly than a realloc
  //but alas c++ is used to be standard
  MemResize(&psp, memsize, ipsp[nnode+gnode]);
  
  checked = new Bool[nnode+gnode];
  
  //init checked list to false
  for(i = 0; i < nnode+gnode; i++){
    checked[i] = 0;
  }

  //perform sanity checks - if node 2 is in node 1's list then node 1 
  //must also be in node 2's list
  for(i = 0; i < nnode+gnode; i++){
    count1 = 0;
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(j = indx1; j < indx2; j++){
      //if reverse map on node isn't fully verified.. check it
      if(checked[psp[j]] == 0){
	indx3 = ipsp[psp[j]];
	indx4 = ipsp[psp[j]+1];
	for(k = indx3; k < indx4; k++){
	  if(psp[k] == i){
	    count1++;
	    break;
	  }
	}
      }
      //if reverse map on node is verified.. we've already checked it
      else{
	count1++;
      }
    }
    if(count1 != (indx2-indx1)){
      std::cout << "MESH UTILITY: WARNING!!! Psp[] map failed sanity check on node " 
		<< i << std::endl;
    }
    else{
      checked[i] = 1;
    }
  }
  
  delete [] checked;
  delete [] cli;
  delete [] cli_d;

  return (0);
}

template <class Type>
Int Mesh<Type>::BuildEl2el()
{
  std::cout << "MESH UTILITY: Building volume element to surrounding element map... " << std::endl;

  el2el = NULL;
  iel2el = NULL;

  Int e, i, j, k, jj;
  Int inode;
  Int gelem, jelem;
  Int indx1, indx2;
  Int size = 1*nelem[TRI] + 1*nelem[QUAD] + 4*nelem[TET]
    + 5*nelem[PYRAMID] + 5*nelem[PRISM] + 6*nelem[HEX];
  Int size2;

  Int knodes [4];
  Int check, face;
  Int facenum [] = {1, 1, 4, 5, 5, 6};
  //nodes contained in faces
  Int fl[] =
    {0, 1, 2, //Tri 
     0, 1, 2, 3, //Quad 
     0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2,         //Tet
     0, 1, 2, 3, 0, 4, 1, 0, 3, 4, 1, 4, 2, 2, 4, 3,   //Pyramid 
     0, 1, 2, 0, 3, 4, 1, 0, 2, 5, 3, 1, 4, 5, 2, 3, 5, 4,       //Prism
     0, 1, 2, 3, 0, 4, 5, 1, 0, 3, 7, 4, 1, 5, 6, 2, 2, 6, 7, 3, 4, 7, 6, 5 //Hex
    };
  Int** fli = NULL;
  Int* fli_d = NULL;

  fli = new Int*[6];
  fli_d = new Int[23];

  el2el = new Int[size];
  iel2el = new Int[lnelem+1];

  if(el2el == NULL || iel2el == NULL || fli == NULL || fli_d == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildEl2el()" << std::endl;
    return (-1);
  }
  //setup CRS indexes by hand for fli
  i = 0;
  fli[TRI] = &fli_d[i];
  i += facenum[TRI];
  fli[QUAD] = &fli_d[i];
  i += facenum[QUAD];
  fli[TET] = &fli_d[i];
  i += facenum[TET];
  fli[PYRAMID] = &fli_d[i];
  i += facenum[PYRAMID];
  fli[PRISM] = &fli_d[i];
  i += facenum[PRISM];
  fli[HEX] = &fli_d[i];

  //fill in CRS indexes
  fli[TRI][0] = 0;
  fli[QUAD][0] = 3;
  fli[TET][0] = 7;
  fli[TET][1] = 10;
  fli[TET][2] = 13;
  fli[TET][3] = 16;
  fli[PYRAMID][0] = 19;
  fli[PYRAMID][1] = 23;
  fli[PYRAMID][2] = 26;
  fli[PYRAMID][3] = 29;
  fli[PYRAMID][4] = 32;
  fli[PRISM][0] = 35;
  fli[PRISM][1] = 38;
  fli[PRISM][2] = 42;
  fli[PRISM][3] = 46;
  fli[PRISM][4] = 50;
  fli[HEX][0] = 53;
  fli[HEX][1] = 57;
  fli[HEX][2] = 61;
  fli[HEX][3] = 65;
  fli[HEX][4] = 69;
  fli[HEX][5] = 73;
  //set last value, not really 7 faces in a hex
  fli[HEX][6] = 77;

  //setup offsets for iel2el
  size2 = 0;
  j = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    e = element.GetType();
    iel2el[j] = size2;
    size2 += facenum[e];
    j++;
  }
  iel2el[lnelem] = size2;
  
  //init el2el map for use in parallel
  //i.e. a negative number equals crossing
  //a non-local connection
  for(i = 0; i < size; i++){
    el2el[i] = -1;
  }

  //load el2el map
  gelem = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes = NULL;
    Int nn = element.GetNodes(&nodes);
    e = element.GetType();
    for(jj = 0; jj < nn; ++jj){
      inode = nodes[jj];
      indx1 = ielsp[inode];
      indx2 = ielsp[inode+1];
      for(j = indx1; j < indx2; j++){
	jelem = elsp[j];
	if(jelem != gelem){
	  Element<Type>& element2 = *elementList[jelem];
	  check = CheckSharedFace(element, element2, knodes);
	  if(check > 0){
	    face = CheckFaceLocation(fli, fl, e, facenum[e], knodes, check);
	    if(face < 0){
	      std::cerr << "WARNING!!! Face connectivity not found -- element " << 
		gelem << " likely broken, type " << e << std::endl;
	    }
	    else{
	      el2el[iel2el[gelem] + face] = jelem;
	    }
	  }
	}
      }
    }
    gelem++;
  }
  Bool* checked = new Bool[lnelem];
  Int indx3, indx4;
  Int count1;
  //init checked list to false
  for(i = 0; i < lnelem; i++){
    checked[i] = 0;
  }
  //perform sanity checks - if elem 2 is in elem 1's list then elem 1 
  //must also be in elem 2's list
  for(i = 0; i < lnelem; i++){
    count1 = 0;
    indx1 = iel2el[i];
    indx2 = iel2el[i+1];
    for(j = indx1; j < indx2; j++){
      //if connection is to another domain set as okay
      //this technically only means the connecting element was never found
      if(el2el[j] == -1){
	count1++;
	continue;
      }
      //if reverse map on elem isn't fully verified.. check it
      if(checked[el2el[j]] == 0){
	indx3 = iel2el[el2el[j]];
	indx4 = iel2el[el2el[j]+1];
	for(k = indx3; k < indx4; k++){
	  if(el2el[k] == i){
	    count1++;
	    break;
	  }
	}
      }
      //if reverse map on elem is verified.. we've already checked it
      else{
	count1++;
      }
    }
    if(count1 != (indx2-indx1)){
      std::cout << "MESH UTILITY: WARNING!!! El2el[] map failed sanity check on element " 
		<< i << std::endl;
    }
    else{
      checked[i] = 1;
    }
  }
  delete [] checked;
  delete [] fli;
  delete [] fli_d;
  
  return (0);
}

template <class Type>
Int Mesh<Type>::BuildSel2el()
{
  std::cout << "MESH UTILITY: Building surface element to surrounding element map... " << std::endl;

  Int size = 3*nelem[TRI] + 4*nelem[QUAD];
  Int gelem, jelem;
  Int inode;
  Int indx1, indx2;
  Int check, edgeid;
  Int knodes[8];
  Int edgenum[] = {3, 4};
  Int el[] = 
    {0, 1, 1, 2, 2, 0, //Tri Edges
     0, 1, 1, 2, 2, 3, 3, 0  //Quad Edges
    };
  Int** eli = NULL;
  Int* elid = NULL;
  isel2el = new Int[lnelem+1];
  sel2el = new Int[size];
  
  eli = new Int*[2];
  elid = new Int[14];
 
  if(sel2el == NULL || isel2el == NULL || eli == NULL || elid == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildSel2el()" << std::endl;
    return (-1);
  }

  //set up edge lookup table
  eli[TRI] = &elid[0];
  eli[QUAD] = &elid[edgenum[TRI]];

  eli[TRI][0] = 0;
  eli[TRI][1] = 2;
  eli[TRI][2] = 4;
  eli[QUAD][0] = 6;
  eli[QUAD][1] = 8;
  eli[QUAD][2] = 10;
  eli[QUAD][3] = 12;
  //set last pointer
  eli[QUAD][4] = 14;

  //setup offsets for isel2el
  Int size2 = 0;
  gelem  = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int e = element.GetType();
    isel2el[gelem] = size2;
    if(e == TRI || e == QUAD){
      size2 += edgenum[e];
    }
    gelem++;
  }
  isel2el[lnelem] = size2;
  
  //init sel2el map for use in parallel
  //i.e. a negative number equals crossing
  //a non-local connection
  for(Int i = 0; i < size; i++){
    sel2el[i] = -1;
  }

  gelem = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int e = element.GetType();
    if(e == TRI || e == QUAD){
      Int* nodes = NULL;
      Int nnodes = element.GetNodes(&nodes);
      for(Int i = 0; i < nnodes; i++){
	inode = nodes[i];
	indx1 = ielsp[inode];
	indx2 = ielsp[inode+1];
	for(Int j = indx1; j < indx2; j++){
	  jelem = elsp[j];
	  Element<Type>& element2 = *elementList[jelem];
	  Int e2 = element2.GetType();
	  //ensure element is a boundary element 
	  if(jelem != gelem && (e2 == TRI || e2 == QUAD)){
	    check = CheckSharedEdge(element, element2, knodes);
	    if(check > 0){
	      edgeid = CheckEdgeLocationSurf(eli, el, e, edgenum[e], knodes, check);
	      if(edgeid < 0){
		std::cerr << "WARNING!!! Edge connectivity not found -- element " << 
		  gelem << " to element " << jelem << " likely broken, type " << e << std::endl;
	      }
	      else{
		sel2el[isel2el[gelem] + edgeid] = jelem;
	      }
	    }
	  }
	}
      }
    }
    gelem++;
  }

  Bool* checked = new Bool[lnelem];
  Int indx3, indx4;
  Int count1;
  //init checked list to false
  for(Int i = 0; i < lnelem; i++){
    checked[i] = 0;
  }
  //perform sanity checks - if elem 2 is in elem 1's list then elem 1 
  //must also be in elem 2's list
  for(Int i = 0; i < lnelem; i++){
    count1 = 0;
    indx1 = isel2el[i];
    indx2 = isel2el[i+1];
    for(Int j = indx1; j < indx2; j++){
      //if connection is to another domain set as okay
      //this technically only means the connecting element was never found
      if(sel2el[j] == -1){
	count1++;
	continue;
      }
      //if reverse map on elem isn't fully verified.. check it
      if(checked[sel2el[j]] == 0){
	indx3 = isel2el[sel2el[j]];
	indx4 = isel2el[sel2el[j]+1];
	for(Int k = indx3; k < indx4; k++){
	  if(sel2el[k] == i){
	    count1++;
	    break;
	  }
	}
      }
      //if reverse map on elem is verified.. we've already checked it
      else{
	count1++;
      }
    }
    if(count1 != (indx2-indx1)){
      std::cout << "MESH UTILITY: WARNING!!! Sel2el[] map failed sanity check on element " 
		<< i << std::endl;
    }
    else{
      checked[i] = 1;
    }
  }

  delete [] checked;
  delete [] eli;
  delete [] elid;

  return (0);
}

template <class Type>
Int Mesh<Type>::BuildEdges()
{
  std::cout << "MESH UTILITY: Building edge map... " << std::endl;
  
  Int i, j, k;
  Int jj;
  Int indx1, indx2;
  Int inode, gelem;
  Int ghostEdgeCount = 0;

  nedge = 0;
  //guess at the number of boundary edges
  nbedge = 3*nelem[TRI] + 4*nelem[QUAD];
  edges = new Edges<Type>[ipsp[nnode]];

  if(edges == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildEdges()" << std::endl;
    return (-1);
  }
  
  //build internal edges
  for(i = 0; i < nnode; i++){
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(j = indx1; j < indx2; j++){
      if(psp[j] >= nnode){
	ghostEdgeCount++;
	continue;
      }
      if(psp[j] > i){
	edges[nedge].n[0] = i;
	edges[nedge].n[1] = psp[j];
	nedge++;
      }
    }
  }

  bedges = new HalfEdges<Type>[nbedge+ghostEdgeCount];
  if(bedges == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildEdges()" << std::endl;
    return (-1);
  }

  //build boundary edges first, this puts bedges first in the list
  //which allows us to directly reference (w/o offset) the edge id for 
  //information stored only on boundaries like yplus, etc.
  k = 0;
  jj = nnode+gnode;
  gelem = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int type = element.GetType();
    if(type == TRI || type == QUAD){
      Int* nodes;
      Int nnodes = element.GetNodes(&nodes);
      for(j = 0; j < nnodes; j++){
	inode = nodes[j];
	//do not create a boundary edge if the boundary face
	//is split and we are on a non-local node
	if(IsInteriorNode(inode)){
	  bedges[k].elem = gelem;
	  //first node is part of the mesh
	  //second node is phantom/ghost node
	  //with no coordinates
	  bedges[k].n[0] = inode;
	  bedges[k].n[1] = jj;
	  bedges[k].factag = element.GetFactag();
	  k++;
	  jj++;
	}
      }
    }
    gelem++;
  }

  //set the true count for the number of BC edges and the number of BC nodes
  nbedge = k;
  nbnode = k;
  //build ghost edges (parallel)
  for(i = 0; i < nnode; i++){
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(j = indx1; j < indx2; j++){
      //this check implies that a ghost node is
      //attached to an interior node...
      if(psp[j] >= nnode){
	bedges[k].elem = -1;
	bedges[k].n[0] = i;
	bedges[k].n[1] = psp[j];
	//note that this is dependent on bc.h
	//though not explicitly done here...
	//be careful changing the factag
	bedges[k].factag = 0;
	k++;
      }
    }
  }

  //set memory init flag
  edgeInit = 1;
  ngedge = ghostEdgeCount;

 //TODO: resize memory to exact number of bedges used

  //allocate cg memory since we have to have an accurate count for nbnode
  //first and this is the first place we can make that statement
 
  //center of gravity -- one per cv * dimensions and one per boundary face edge
  cg = new Type[(nnode+gnode+nbnode)*3];

  return (0);
}

template <class Type>
Int Mesh<Type>::BuildEsp()
{
  Int i, j, ii;
  Int indx1, indx2;
  std::cout << "MESH UTILITY: Building point to edge map... " << std::endl;
  
  esp = NULL;
  iesp = NULL;
  besp = NULL;
  ibesp = NULL;

  iesp = new Int[nnode+gnode+1];
  ibesp = new Int[nnode+gnode+1];

  if(iesp == NULL || ibesp == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildEsp()" << std::endl;
    return (-1);
  }
  
  //init list
  for(i = 0; i <= nnode+gnode; i++){
    ibesp[i] = 0;
    iesp[i] = 0;
  }
  for(i = 0; i < nnode+gnode; i++){
    indx1 = iselsp[i];
    indx2 = iselsp[i+1];
    //add one bedge for each surface neighbor element
    ibesp[i+1] += (indx2-indx1);
  }
  for(i = 0; i < nnode+gnode; i++){
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(j = indx1; j < indx2; j++){
      //if connected to a ghost node add a bedge
      if(psp[j] >= nnode){
	ibesp[i+1]++;
      }
    }
  }
  //setup offsets
  for(i = 1; i <= nnode+gnode; i++){
    ibesp[i] += ibesp[i-1];
  } 
  besp = new Int[ibesp[nnode+gnode]];
 
  if(besp == NULL){
    std::cerr << "MESH UTILITY: Memory allocation failed Mesh::BuildEsp()" << std::endl;
    return (-1);
  }

  for(i = 0; i < nbedge+ngedge; i++){
    ii = bedges[i].n[0];
    besp[ibesp[ii]] = i;
    ibesp[ii]++;
  }
  //shift indices in CRS back
  for(i = nnode+gnode; i > 0; i--){
    ibesp[i] = ibesp[i-1];
  }
  ibesp[0] = 0;


  //calculate offsets for edge surrounding point
  for(i = 0; i < nnode+gnode; i++){
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(j = indx1; j < indx2; j++){
      //if point connected is not a ghost node
      if(psp[j] < nnode){
	iesp[i+1]++;
      }
    }
  }
  for(i = 1; i <= nnode+gnode; i++){
    iesp[i] += iesp[i-1];
  } 
  iesp[0] = 0;
  
  esp = new Int[iesp[nnode+gnode]];
  for(i = 0; i < nedge; i++){
    esp[iesp[edges[i].n[0]]] = i;
    iesp[edges[i].n[0]]++;
    esp[iesp[edges[i].n[1]]] = i;
    iesp[edges[i].n[1]]++;
  }

  //shift CRS indexes back
  for(i = nnode+gnode; i > 0; i--){
    iesp[i] = iesp[i-1];
  }
  iesp[0] = 0;
  
  return (0);
}

template <class Type>
Bool Mesh<Type>::IsElemParallel(Int gelem)
{
  Element<Type>& element = *elementList[gelem];
  Int* nodes = NULL;
  Int nnodes = element.GetNodes(&nodes);

  for(Int i = 0; i < nnodes; i++){
    if(IsGhostNode(nodes[i])){
      return true;
    }
  }
  return false;
}

template <class Type>
void Mesh<Type>::GetParallelTuple(Int node, Int* owningProcess, Int* localId)
{
  if(!IsGhostNode(node)){
    *owningProcess = -1;
    *localId = -1;
    return;
  }
  else{
    node -= nnode;
    *localId = gNodeLocalId[node];
    *owningProcess = gNodeOwner[node];
    return;
  }
}



/******************************************/
//Find repeats in a list, reorder list with
//repeats at the beginning and return
//number of unique values
/******************************************/
template <class Type>
Int Mesh<Type>::EliminateRepeats(Int* array, Int n)
{
  Int i,j;
  Int unique;
  //assume the local list is small and won't 
  //overflow the stack -- likely safe here
  Bool flag[n]; 

  for(i = 0; i < n; i++){
    flag[i] = 1;
  }
  //mark values which are not unique
  //with boolean false
  for(i = 0; i < n; i++){
    //make sure we haven't marked value i
    //for deletion -- avoid repetitive checks
    if(flag[i] == 1){
      for(j = i+1; j < n; j++){
	if(array[i] == array[j]){
	  flag[j] = 0;
	}
      }
    }
  }
  //move unique values to beginning of list
  //first value always checks as unique
  unique = 1;
  for(i = 1; i < n; i++){
    if(flag[i] == 1){
      array[unique] = array[i];
      unique++;
    }
  }
  return (unique);
}

/********************************************/
//Checks for shared face between two elements
//via a node list compare.. returns faces if 
//shares a face (3+ nodes) and 0 if they do not
//This is expensive -- room for improvement!!
//knodes[] is set to the integer value 0-4 of
//nodes that matched
/********************************************/
template <class Type>
Int Mesh<Type>::CheckSharedFace(Element<Type>& elem1, Element<Type>& elem2, Int* knodes)
{
  Int nodes = 0;
  Int* nodes1 = NULL;
  Int* nodes2 = NULL;
  Int nn1 = elem1.GetNodes(&nodes1);
  Int nn2 = elem2.GetNodes(&nodes2);
  
  for(Int i = 0; i < nn1; i++){
    for(Int j = 0; j < nn2; j++){
      if(nodes1[i] == nodes2[j]){
	knodes[nodes] = i;
	nodes++;
      }
    }
  }
  if(nodes >= 3){
    return (nodes);
  }
  else{
    return (0);
  }
}

/********************************************/
//Same function as above but only checks for
//a shared edge -- returns the points of the
//edge in knodes
/********************************************/
template <class Type>
Int Mesh<Type>::CheckSharedEdge(Element<Type>& elem1, Element<Type>& elem2, Int* knodes)
{
  Int nodes = 0;
  Int* nodes1 = NULL;
  Int* nodes2 = NULL;
  Int nn1 = elem1.GetNodes(&nodes1);
  Int nn2 = elem2.GetNodes(&nodes2);
  
  for(Int i = 0; i < nn1; i++){
    for(Int j = 0; j < nn2; j++){
      if(nodes1[i] == nodes2[j]){
	knodes[nodes] = i;
	nodes++;
	break;
      }
    }
  }
  if(nodes >= 2){
    return (nodes);
  }
  else{
    return (0);
  }
}

/******************************************/
//given lookup tables fli and fl, returns
//the face number for a list of nodes
//given in arguments
/******************************************/
template <class Type>
Int Mesh<Type>::CheckFaceLocation(Int** fli, Int* fl, Int etype, Int facenum, Int* knodes, Int nodes)
{
  for(Int i = 0; i < facenum; i++){
    Int count = 0;
    Int indx1 = fli[etype][i];
    Int indx2 = fli[etype][i+1];
    for(Int j = indx1; j < indx2; j++){
      for(Int k = 0; k < nodes; k++){
	if(knodes[k] == fl[j]){
	  count++;
	}
	if(count >= 3){
	  return(i);
	}
      }
    }
  }
  return(-1);
}

/***************************************/
//given lookup table eli and el, returns
//the edge number for a list of nodes 
//given. Only applicable to 2D, in 3D
//more than one edge could be shared and
//a list should be returned
/***************************************/
template <class Type>
Int Mesh<Type>::CheckEdgeLocationSurf(Int** eli, Int* el, Int etype, Int ednum, Int* knodes, Int nodes)
{
  for(Int i = 0; i < ednum; i++){
    Int count = 0;
    Int indx1 = eli[etype][i];
    Int indx2 = eli[etype][i+1];
    for(Int j = indx1; j < indx2; j++){
      for(Int k = 0; k < nodes; k++){
	if(knodes[k] == el[j]){
	  count++;
	}
	if(count >= 2){
	  return(i);
	}
      }
    }
  }
  return(-1);
}

/*****************************************/
//Given a particular edges nodes and element (g)
//check if that edge is part of given
//element, return 1 if true, 0 if false
/*****************************************/
template <class Type>
Int Mesh<Type>::EdgeInElem(Int n1, Int n2, Int jelem)
{
  Int count = 0;
  Int* nodes = NULL;
  Int ennode = elementList[jelem]->GetNodes(&nodes);

  for(Int i = 0; i < ennode; ++i){
    Int jnode = nodes[i];
    if(jnode == n1 || jnode == n2){
      count++;
    }
  }
  if(count == 2){
    return (1);
  }
  return (0);
}


/*****************************************/
//Function which calls necessary functions
//to calculate grid metrics
/*****************************************/
template <class Type>
Int Mesh<Type>::CalcMetrics()
{
  Int err = 0;
  std::cout << "MESH UTILITY: Calculating geometry metrics" << std::endl;
  err = CalcExtents();
  if(err != 0){
    return(err);
  }
  err = CalcAreasVolumes();
  if(err != 0){
    return (err);
  }
  metricsCalc = 1;
  err = CheckForClosure();
  if(err != 0){
    return (err);
  }
  return (0);
}

template <class Type>
Int Mesh<Type>::CalcExtents()
{
  Int i;
  Real maxX, maxY, maxZ;
  Real minX, minY, minZ;

  MPI_Datatype mpit;

  //silence warning
  maxX = maxY = maxZ = minX = minY = minZ = 0.0;

  //get mpi datatype to send
  Real tmp = 0.0;
  mpit = MPI_GetType(tmp);

  minX = maxX = real(xyz[0*3 + 0]);
  minY = maxY = real(xyz[0*3 + 1]);
  minZ = maxZ = real(xyz[0*3 + 2]);
  for(i = 1; i < nnode+gnode; i++){
    maxX = real(MAX(maxX, xyz[i*3 + 0]));
    maxY = real(MAX(maxY, xyz[i*3 + 1]));
    maxZ = real(MAX(maxZ, xyz[i*3 + 2]));
    minX = real(MIN(minX, xyz[i*3 + 0]));
    minY = real(MIN(minY, xyz[i*3 + 1]));
    minZ = real(MIN(minZ, xyz[i*3 + 2]));
  }

  Real* realExMax = (Real*)alloca(sizeof(Real)*3);
  Real* realExMin = (Real*)alloca(sizeof(Real)*3);

  realExMax[0] = maxX;
  realExMax[1] = maxY;
  realExMax[2] = maxZ;
  realExMin[0] = minX;
  realExMax[1] = minY;
  realExMax[2] = minZ;

  MPI_Allreduce(MPI_IN_PLACE, realExMax, 3, mpit, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, realExMin, 3, mpit, MPI_MIN, MPI_COMM_WORLD);

  extentsMax[0] = realExMax[0];
  extentsMax[1] = realExMax[1];
  extentsMax[2] = realExMax[2];
  extentsMin[0] = realExMin[0];
  extentsMin[1] = realExMin[1];
  extentsMin[2] = realExMin[2];

  std::cout << "MESH UTILITY: Maximum mesh extents: " << extentsMax[0] << " " 
	    << extentsMax[1] << " " << extentsMax[2] <<std::endl;
  std::cout << "MESH UTILITY: Minimum mesh extents: " << extentsMin[0] << " " 
	    << extentsMin[1] << " " << extentsMin[2] <<std::endl;

  return(0);
}

template <class Type>
Int Mesh<Type>::CalcAreasVolumes()
{
  Int i, j, k, jj, kk, nn;
  Int indx1, indx2, indx3, indx4;
  Int n1, n2;
  Int size = 12;
  Int count;
  Int* elems = new Int[20];
  Int nloc[2];
  Int jelem, e;
  Int jnode;
  Int hits;
  Int orientation = -1;
  Int eid;
  Int face1 = -1;
  Int face2 = -1;
  Type edgept[3];
  Type face1pt[3];
  Type face2pt[3];
  Type elempt[3];
  Int* nodes;
  Int nodef[4];
  Type area[3];
  //temp vol
  Type tv[4];
  //temp cg
  Type tcg[4*3];
  Type meshCg[3];
  Type sum = 0.0;
 
  //number of faces in each element type
  Int facenum [] = {1, 1, 4, 5, 5, 6};
  //don't consider boundary elements at this time
  Int numedges[] = {0, 0, 6, 8, 9, 12};
  //edge ordering -- two values consecutively represent
  //an edge and its preferred orientation
  //i.e. first edge in a tet (eid = 0) is oriented 0,1 (1st 2 ints)
  Int eorder[] =
    {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, //Tet
     0, 1, 0, 3, 0, 4, 1, 2, 1, 4, 2, 3, 2, 4, 3, 4, //Pyramid
     0, 1, 0, 2, 0, 3, 1, 2, 1, 4, 2, 5, 3, 4, 3, 5, 4, 5, //Prism
     0, 1, 0, 3, 0, 4, 1, 2, 1, 5, 2, 3, 2, 6, 3, 7, 4, 5, 4, 7, 5, 6, 6, 7//Hex      
    };
  Int** ieorder = new Int*[6];
  Int* ieorderd = new Int[36];

  i = 0;
  ieorder[TRI] = NULL;
  ieorder[QUAD] = NULL;
  ieorder[TET] = &ieorderd[i];
  i += numedges[TET];
  ieorder[PYRAMID] = &ieorderd[i];
  i += numedges[PYRAMID];
  ieorder[PRISM] = &ieorderd[i];
  i += numedges[PRISM];
  ieorder[HEX] = &ieorderd[i];

  i = 0;
  for(e = TET; e <= HEX; e++){
    for(j = 0; j < numedges[e]; j++){
      ieorder[e][j] = i;
      i += 2;
    }
  }
  ieorder[HEX][numedges[HEX]] = i;

  //ordered face to edge such that the first face
  //list is to the left of the oriented edge when
  //looking from the inside of the element
  //i.e. faces l,r of edge 0 in a Tet are 0,1 (1st 2 ints)
  Int e2f[] =
    {0, 1, 2, 0, 1, 2, 0, 3, 3, 1, 2, 3, //Tet
     0, 1, 2, 0, 1, 2, 0, 3, 3, 1, 0, 4, 4, 3, 2, 4, //Pyramid
     0, 1, 2, 0, 1, 2, 0, 3, 3, 1, 2, 3, 1, 4, 4, 2, 3, 4, //Prism
     0, 1, 2, 0, 1, 2, 0, 3, 3, 1, 0, 4, 4, 3, 2, 4, 1, 5, 5, 2, 3, 5, 4, 5 //Hex     
    };
  
  Int** ie2f = new Int*[6];
  Int* ie2fd = new Int[36];
  
  i = 0;
  ie2f[TRI] = NULL;
  ie2f[QUAD] = NULL;
  ie2f[TET] = &ie2fd[i];
  i += numedges[TET];
  ie2f[PYRAMID] = &ie2fd[i];
  i += numedges[PYRAMID];
  ie2f[PRISM] = &ie2fd[i];
  i += numedges[PRISM];
  ie2f[HEX] = &ie2fd[i];

  i = 0;
  for(e = TET; e <= HEX; e++){
    for(j = 0; j < numedges[e]; j++){
      ie2f[e][j] = i;
      i += 2;
    }
  }
  ie2f[HEX][numedges[HEX]] = i;
 
  //nodes contained in faces
  Int fl[] =
    {0, 1, 2, //Tri 
     0, 1, 2, 3, //Quad 
     0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2,         //Tet
     0, 1, 2, 3, 0, 4, 1, 0, 3, 4, 1, 4, 2, 2, 4, 3,   //Pyramid 
     0, 1, 2, 0, 3, 4, 1, 0, 2, 5, 3, 1, 4, 5, 2, 3, 5, 4,       //Prism
     0, 1, 2, 3, 0, 4, 5, 1, 0, 3, 7, 4, 1, 5, 6, 2, 2, 6, 7, 3, 4, 7, 6, 5 //Hex
    };
  Int** fli = NULL;
  Int* fli_d = NULL;
  fli = new Int*[6];
  fli_d = new Int[23];
  //setup CRS indexes by hand for fli
  i = 0;
  fli[TRI] = &fli_d[i];
  i += facenum[TRI];
  fli[QUAD] = &fli_d[i];
  i += facenum[QUAD];
  fli[TET] = &fli_d[i];
  i += facenum[TET];
  fli[PYRAMID] = &fli_d[i];
  i += facenum[PYRAMID];
  fli[PRISM] = &fli_d[i];
  i += facenum[PRISM];
  fli[HEX] = &fli_d[i];
  //fill in CRS indexes
  fli[TRI][0] = 0;
  fli[QUAD][0] = 3;
  fli[TET][0] = 7;
  fli[TET][1] = 10;
  fli[TET][2] = 13;
  fli[TET][3] = 16;
  fli[PYRAMID][0] = 19;
  fli[PYRAMID][1] = 23;
  fli[PYRAMID][2] = 26;
  fli[PYRAMID][3] = 29;
  fli[PYRAMID][4] = 32;
  fli[PRISM][0] = 35;
  fli[PRISM][1] = 38;
  fli[PRISM][2] = 42;
  fli[PRISM][3] = 46;
  fli[PRISM][4] = 50;
  fli[HEX][0] = 53;
  fli[HEX][1] = 57;
  fli[HEX][2] = 61;
  fli[HEX][3] = 65;
  fli[HEX][4] = 69;
  fli[HEX][5] = 73;
  //set last value, not really 7 faces in a hex
  fli[HEX][6] = 77;

  //zero volumes/areas in mesh
  for(i = 0; i < nnode; i++){
    vol[i] = 0.0;
    cg[i*3 + 0] = cg[i*3 + 1] = cg[i*3 + 2] = 0.0;
  }
  for(i = 0; i < nedge; i++){
    edges[i].a[0] = edges[i].a[1] = edges[i].a[2] = edges[i].a[3] = 0.0;
  }
  for(i = 0; i < ngedge+nbedge; i++){
    bedges[i].a[0] = bedges[i].a[1] = bedges[i].a[2] = bedges[i].a[3] = 0.0;
  }

  /***************************************/
  //Close off portions of CVs with edges
  /***************************************/
  for(i = 0; i < nedge; i++){
    count = 0;
    n1 = edges[i].n[0];
    n2 = edges[i].n[1];
    indx1 = ielsp[n2];
    indx2 = ielsp[n2+1];
    for(k = indx1; k < indx2; k++){
      jelem = elsp[k];
      Element<Type>& element = *elementList[jelem];
      Int jtype = element.GetType();
      //do not include bfaces -- redundant and unneccessary
      //since a bface lies coincident with some volume
      //element face anyway
      if(!(jtype == TRI || jtype == QUAD)){
	if(EdgeInElem(n1, n2, jelem)){
	  elems[count] = jelem;
	  count++;
	  if(count == size){
	    MemResize(&elems, size, 2*size);
	    size += size;
	  }
	}
      }
    }
    for(j = 0; j < count; j++){
      hits = 0;
      jelem = elems[j];
      Element<Type>& element = *elementList[jelem];
      Int nnodes = element.GetNodes(&nodes);
      Int jtype = element.GetType();
      for(k = 0; k < nnodes; k++){
	jnode = nodes[k];
	if(n1 == jnode){
	  nloc[0] = k;
	  hits++;
	}
	else if(n2 == jnode){
	  nloc[1] = k;
	  hits++;
	}
	if(hits == 2){
	  orientation = CheckEdgeOrientation(ieorder, eorder, jtype, numedges[jtype], nloc, eid);
	  if(orientation == 0){
	    std::cerr << "WARNING: cannot determine edge orientation -- Mesh::CalcMetrics() " 
		      << "for element " << jelem << std::endl;
	  }
	  else{
	    face1 = e2f[ie2f[jtype][eid]];
	    face2 = e2f[ie2f[jtype][eid]+1];
	  }
	  break;
	}
      }
      if(hits != 2){
	std::cerr << "WARNING: CalcMetrics() -- edge orientation not found, element " << jelem << std::endl;
      }

      //calculate edge midpoint
      Centroid(&xyz[n1*3 + 0], &xyz[n2*3 + 0], edgept);
      
      //calculate face1 centroid
      indx1 = fli[jtype][face1];
      indx2 = fli[jtype][face1+1];
      jj = 0;
      for(kk = indx1; kk < indx2; kk++){
	nodef[jj] = nodes[fl[kk]];
	jj++;
      }
      if(indx2 - indx1 == 3){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 face1pt);
      }
      else if(indx2 - indx1 == 4){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 &xyz[nodef[3]*3], face1pt);
      }
      else{
	std::cerr << "WARNING!!! CalcMetrics() -- odd face detected with " 
		  << indx2-indx1 << " nodes - in element " << jelem << std::endl;
      }
      
      //calculate face2 centroid
      indx3 = fli[jtype][face2];
      indx4 = fli[jtype][face2+1];
      jj = 0;
      for(kk = indx3; kk < indx4; kk++){
	nodef[jj] = nodes[fl[kk]];
	jj++;
      }
      if(indx4 - indx3 == 3){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 face2pt);
      }
      else if(indx4 - indx3 == 4){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 &xyz[nodef[3]*3], face2pt);
      }
      else{
	std::cerr << "WARNING!!! CalcMetrics() -- odd face detected with " 
		  << indx4-indx3 << " nodes in element " << jelem << std::endl;
      }

      //calculate element centroid
      ElementCentroid(nodes, xyz, jtype, elempt);

      //calculate area of partial face
      CalcTriArea(elempt, face2pt, edgept, area);
      Scale((Type)orientation, area, 3);
      Add(edges[i].a, area, edges[i].a);
 
      CalcTriArea(edgept, face1pt, elempt, area);
      Scale((Type)orientation, area, 3);
      Add(edges[i].a, area, edges[i].a);

      //calculate volume of partial CV
      tv[0] = (Type)orientation * CalcTetVol(&xyz[n1*3], elempt, face2pt, edgept);
      tv[1] = (Type)orientation * CalcTetVol(&xyz[n1*3], face1pt, elempt, edgept);
      tv[2] = (Type)orientation * CalcTetVol(&xyz[n2*3], elempt, face1pt, edgept);
      tv[3] = (Type)orientation * CalcTetVol(&xyz[n2*3], face2pt, elempt, edgept);

      //calculate cg of each partial CV
      Centroid(&xyz[n1*3], elempt, face2pt, edgept, &tcg[0*3]);
      Centroid(&xyz[n1*3], face1pt, elempt, edgept, &tcg[1*3]);
      Centroid(&xyz[n2*3], elempt, face1pt, edgept, &tcg[2*3]);
      Centroid(&xyz[n2*3], face2pt, elempt, edgept, &tcg[3*3]);

      //weight each partial cg coordinate by volume of partial CV
      for(nn = 0; nn < 4; nn++){
	for(kk = 0; kk < 3; kk++){
	  tcg[nn*3 + kk] *= tv[nn];
	}
      }
      for(kk = 0; kk < 4; kk++){
	if(real(tv[kk]) < 0.0 ){//&& Abs(tv[kk]) > 1.0e-15){
	  std::cout << "Negative volume detected in element " << jelem << " of type " << jtype 
		    << " -- volume " << tv[kk] << " # " << kk << std::endl;
	}
      }
      vol[n1] += tv[0];
      vol[n1] += tv[1];
      vol[n2] += tv[2];
      vol[n2] += tv[3];
      for(kk = 0; kk < 3; kk++){
	cg[n1*3 + kk] += tcg[0*3 + kk];
	cg[n1*3 + kk] += tcg[1*3 + kk];
      }
      for(kk = 0; kk < 3; kk++){
	cg[n2*3 + kk] += tcg[2*3 + kk];
	cg[n2*3 + kk] += tcg[3*3 + kk];
      }
    }
  }

  /***************************************/
  //Close off portions of CVs with gedges
  /***************************************/
  for(i = nbedge; i < ngedge+nbedge; i++){
    count = 0;
    n1 = bedges[i].n[0]; //node in mesh
    n2 = bedges[i].n[1]; //node ghosted in mesh
    indx1 = ielsp[n2];
    indx2 = ielsp[n2+1];
    for(k = indx1; k < indx2; k++){
      jelem = elsp[k];
      Element<Type>& element = *elementList[jelem];
      Int jtype = element.GetType();
      //do not include bfaces -- redundant and unneccessary
      //since a bface lies coincident with some volume
      //element face anyway
      if(!(jtype == TRI || jtype == QUAD)){
	if(EdgeInElem(n1, n2, jelem)){
	  elems[count] = jelem;
	  count++;
	  if(count == size){
	    MemResize(&elems, size, 2*size);
	    size += size;
	  }
	}
      }
    }
    for(j = 0; j < count; j++){
      hits = 0;
      jelem = elems[j];
      Element<Type>& element = *elementList[jelem];
      Int nnodes = element.GetNodes(&nodes);
      Int jtype = element.GetType();
      for(k = 0; k < nnodes; k++){
	jnode = nodes[k];
	if(n1 == jnode){
	  nloc[0] = k;
	  hits++;
	}
	else if(n2 == jnode){
	  nloc[1] = k;
	  hits++;
	}
	if(hits == 2){
	  orientation = CheckEdgeOrientation(ieorder, eorder, jtype, numedges[jtype], nloc, eid);
	  if(orientation == 0){
	    std::cerr << "WARNING: cannot determine edge orientation -- Mesh::CalcMetrics() " 
		      << "for element " << jelem << std::endl;
	  }
	  else{
	    face1 = e2f[ie2f[jtype][eid]];
	    face2 = e2f[ie2f[jtype][eid]+1];
	  }
	  break;
	}
      }
      if(hits != 2){
	std::cerr << "WARNING: CalcMetrics() -- edge orientation not found, element " << jelem << std::endl;
      }

      //calculate edge midpoint
      Centroid(&xyz[n1*3 + 0], &xyz[n2*3 + 0], edgept);
      
      //calculate face1 centroid
      indx1 = fli[jtype][face1];
      indx2 = fli[jtype][face1+1];
      jj = 0;
      for(kk = indx1; kk < indx2; kk++){
	nodef[jj] = nodes[fl[kk]];
	jj++;
      }
      if(indx2 - indx1 == 3){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 face1pt);
      }
      else if(indx2 - indx1 == 4){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 &xyz[nodef[3]*3], face1pt);
      }
      else{
	std::cerr << "WARNING!!! CalcMetrics() -- odd face detected with " 
		  << indx2-indx1 << " nodes - in element " << jelem << std::endl;
      }
      
      //calculate face2 centroid
      indx3 = fli[jtype][face2];
      indx4 = fli[jtype][face2+1];
      jj = 0;
      for(kk = indx3; kk < indx4; kk++){
	nodef[jj] = nodes[fl[kk]];
	jj++;
      }
      if(indx4 - indx3 == 3){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 face2pt);
      }
      else if(indx4 - indx3 == 4){
	Centroid(&xyz[nodef[0]*3], &xyz[nodef[1]*3], &xyz[nodef[2]*3], 
		 &xyz[nodef[3]*3], face2pt);
      }
      else{
	std::cerr << "WARNING!!! CalcMetrics() -- odd face detected with " 
		  << indx4-indx3 << " nodes in element " << jelem << std::endl;
      }

      //calculate element centroid
      ElementCentroid(nodes, xyz, jtype, elempt);

      //calculate area of partial face
      CalcTriArea(elempt, face2pt, edgept, area);
      Scale((Type)orientation, area, 3);
      Add(bedges[i].a, area, bedges[i].a);
 
      CalcTriArea(edgept, face1pt, elempt, area);
      Scale((Type)orientation, area, 3);
      Add(bedges[i].a, area, bedges[i].a);

      //calculate volume of partial CV
      tv[0] = (Type)orientation * CalcTetVol(&xyz[n1*3], elempt, face2pt, edgept);
      tv[1] = (Type)orientation * CalcTetVol(&xyz[n1*3], face1pt, elempt, edgept);

      //calculate cg of each partial CV
      Centroid(&xyz[n1*3], elempt, face2pt, edgept, &tcg[0*3]);
      Centroid(&xyz[n1*3], face1pt, elempt, edgept, &tcg[1*3]);

      //weight each partial cg coordinate by volume of partial CV
      for(nn = 0; nn < 2; nn++){
	for(kk = 0; kk < 3; kk++){
	  tcg[nn*3 + kk] *= tv[nn];
	}
      }
      for(kk = 0; kk < 2; kk++){
	if(real(tv[kk]) < 0.0 ){//&& Abs(tv[kk]) > 1.0e-15){
	  std::cout << "Negative volume detected in element " << jelem << " of type " << jtype 
		    << " -- volume " << tv[kk] << " # " << kk << std::endl;
	}
      }
      vol[n1] += tv[0];
      vol[n1] += tv[1];
      for(kk = 0; kk < 3; kk++){
	cg[n1*3 + kk] += tcg[0*3 + kk];
	cg[n1*3 + kk] += tcg[1*3 + kk];
      }
    }
  }

  /*************************************/
  //Finish closing CVs with bedges
  /************************************/
  //CRS storage for element point neighbor lookups 86 values
  //nodes connected to a particular point in a node
  Int nid = -1;
  Int altnode;
  Int cl[] = 
    {1, 2, 2, 0, 0, 1, //Tri 
     1, 3, 2, 0, 3, 1, 0, 2  //Quad 
    };
  Int** cli = new Int*[2];
  Int* clid = new Int[8];
  //setup lookup table by hand
  i = 0;
  cli[TRI] = &clid[i];
  i += mnode[TRI];
  cli[QUAD] = &clid[i];
  //set CRS indexes
  cli[TRI][0] = 0; 
  cli[TRI][1] = 2;
  cli[TRI][2] = 4;
  cli[QUAD][0] = 6;
  cli[QUAD][1] = 8;
  cli[QUAD][2] = 10;
  cli[QUAD][3] = 12;
  cli[QUAD][4] = 14;

  //loop over only the ghost edges (boundary)
  for(i = 0; i < nbedge; i++){
    n1 = bedges[i].n[0];
    n2 = bedges[i].n[1];
    jelem = bedges[i].elem;
    Element<Type>& element = *elementList[jelem];
    Int nnodes = element.GetNodes(&nodes);
    Int jtype = element.GetType();
    for(k = 0; k < nnodes; k++){
      if(nodes[k] == n1){
	nid = k;
      }
    }
    ElementCentroid(nodes, xyz, jtype, elempt);

    altnode = cl[cli[jtype][nid]];
    altnode = nodes[altnode];
    Centroid(&xyz[n1*3 + 0], &xyz[altnode*3 + 0], edgept);
    CalcTriArea(elempt, &xyz[n1*3], edgept, area);
    Add(&bedges[i].a[0], area, &bedges[i].a[0]);
    //compute cg of first triangle
    Centroid(elempt, &xyz[n1*3], edgept, &tcg[1*3]);
    //weight each partial cg by the area of bedge associated with it
    tv[0] = Magnitude(area);
    for(kk = 0; kk < 3; kk++){
      tcg[0*3 + kk] = tcg[1*3 + kk]*Magnitude(area);
    }
    
    altnode = cl[cli[jtype][nid]+1];
    altnode = nodes[altnode];
    Centroid(&xyz[n1*3 + 0], &xyz[altnode*3 + 0], edgept);
    CalcTriArea(edgept, &xyz[n1*3], elempt, area);
    Add(&bedges[i].a[0], area, &bedges[i].a[0]);
    //compute cg of second triangle
    Centroid(edgept, &xyz[n1*3], elempt, &tcg[2*3]);
    //weight each partial cg by the area of bedge associated with it
    tv[0] += Magnitude(area);
    for(kk = 0; kk < 3; kk++){
      tcg[0*3 + kk] += tcg[2*3 + kk]*Magnitude(area);
    }
    
    //now divide through by the total area to get the weighted cg for 
    //each boundary edge
    for(kk = 0; kk < 3; kk++){
      cg[n2*3 + kk] = tcg[0*3 + kk]/tv[0];
    }
  }

  //cg[] contains the volume weighted coordinates summed for each node
  //divide each by total volume to get actual cg for the CV
  meshCg[0] = meshCg[1] = meshCg[2] = 0.0;
  for(i = 0; i < nnode; i++){
    for(j = 0; j < 3; j++){
      meshCg[j] += cg[i*3 + j];
      cg[i*3 + j] /= vol[i];
    }
  }
  //do a parallel sync for cg's at this point
  p->UpdateGeneralVectors(cg, 3);

  for(i = 0; i < nnode; i++){
    sum += vol[i];
  }
  for(i = 0; i < 3; i++){
    meshCg[i] /= sum;
  }
  std::cout << "MESH UTILITY: Enclosed volume center of mass: " 
	    << meshCg[0] << " " << meshCg[1] << " " << meshCg[2] << std::endl;

  std::cout << "MESH UTILITY: Total enclosed volume: " << sum << std::endl;

  area[0] = area[1] = area[2] = 0.0;
  for(i = 0; i < ngedge+nbedge; i++){
    area[0] += bedges[i].a[0];
    area[1] += bedges[i].a[1];
    area[2] += bedges[i].a[2];
  }
  sum = Magnitude(area);
  std::cout << "MESH UTILITY: Component surface area leakage: " << area[0] 
	    << " " << area[1] << " " << area[2] << std::endl;
  std::cout << "MESH UTILITY: Total surface area leakage: " << sum << std::endl;

  delete [] elems;
  delete [] ieorder;
  delete [] ieorderd;
  delete [] fli;
  delete [] fli_d;
  delete [] ie2f;
  delete [] ie2fd;
  delete [] cli;
  delete [] clid;

  return (0);
}
/****************************************************/
//Given two nodes of an element (locations in element)
//return -1 if in reverse order and 1 if in same order
//also set edgeid to edge number in element --
//needed to find mated faces
//must pass lookup tables as well
//returns 0 on failure....
/****************************************************/
template <class Type>
Int Mesh<Type>::CheckEdgeOrientation(Int** ieorder, Int* eorder, Int etype, Int numedges, 
			       Int* tested, Int& edgeid){
  Int j;
  for(j = 0; j < numedges; j++){
    if(eorder[ieorder[etype][j]] == tested[0] && eorder[ieorder[etype][j]+1] == tested[1]){
      edgeid = j;
      return(1);
    }
    if(eorder[ieorder[etype][j]+1] == tested[0] && eorder[ieorder[etype][j]] == tested[1]){
      edgeid = j;
      return(-1);
    }
  }
  return (0);
}

/**************************************************/
//checks for closure of CVs, returns (-1) if not
//closed and 0 if closed
/**************************************************/
template <class Type>
Int Mesh<Type>::CheckForClosure()
{
  Int i;
  Int n1, n2;
  Bool fail = 0;
  Type smallnum = 1.0e-10;
  Int count = 0;

  Type* areasum = new Type[(nnode)*3];
  for(i = 0; i < nnode; i++){
    areasum[i*3 + 0] = 0.0;
    areasum[i*3 + 1] = 0.0;
    areasum[i*3 + 2] = 0.0;
  }
  for(i = 0; i < nedge; i++){
    n1 = edges[i].n[0];
    n2 = edges[i].n[1];
    Subtract(&areasum[n1*3], edges[i].a, &areasum[n1*3]);
    Add(&areasum[n2*3], edges[i].a, &areasum[n2*3]);
  }
  for(i = 0; i < ngedge+nbedge; i++){
    n1 = bedges[i].n[0];
    Subtract(&areasum[n1*3], bedges[i].a, &areasum[n1*3]);
  }   
  for(i = 0; i < nnode; i++){
    if(real(Magnitude(&areasum[i*3])) > real(smallnum)){
      fail = 1;
      std::cout << "WARNING!!! CV based on node " << i << " not closed -- "
      		<< areasum[i*3 + 0] << " " << areasum[i*3 + 1] << " " 
      		<< areasum[i*3 + 2]<< std::endl;
      count++;
    }
  }
  std::cout << "MESH UTILITY: " << count << " control volumes not watertight... " << std::endl;
  
  delete [] areasum;
  
  //once we have checked and confirmed the closure of the control volumes
  //strip out the magnitude of each and store it in the 4th entry of the
  //edge area vector while normalizing the rest
  for(i = 0; i < nedge; i++){
    edges[i].a[3] = Normalize(edges[i].a, edges[i].a);
  }
  for(i = 0; i < ngedge+nbedge; i++){
    bedges[i].a[3] = Normalize(bedges[i].a, bedges[i].a);
  }

  if(fail == false){
    return (0);
  }
  else{
    return (1);
  }
}

template <class Type>
Int Mesh<Type>::FixBoundaryWindings()
{
  Int k;
  Int indx1, indx2;
  Int jelem;
  Type ne[3], nb[3];
  Int* nodesi;
  Int* nodesj;
  Int nodesitemp[4];
  Type elempt[3], bpt[3];
  Type dot;
  Int count = 0;

  Int gelem = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& elementi = **it;
    Int itype = elementi.GetType();
    if(itype == TRI || itype == QUAD){
      indx1 = iel2el[gelem];
      indx2 = iel2el[gelem+1];
      if(indx2 - indx1 != 1){
	std::cerr << "WARNING!!! Boundary face # " << gelem << " of type " 
		  << itype << " connected to " 
		  << "multiple volume elements Mesh::FixBoundaryWindings()" 
		  << std::endl;
	std::cerr << "Element list: " << std::endl;
	for(k = indx1; k < indx2; k++){
	  std::cerr << el2el[k] << " ";
	}
	std::cerr << std::endl;
	return (-1);
      }
      elementi.GetNodes(&nodesi);
      ElementCentroid(nodesi, xyz, itype, bpt);
      jelem = el2el[indx1];
      Element<Type>& elementj = *elementList[jelem];
      elementj.GetNodes(&nodesj);
      Int jtype = elementj.GetType();
      ElementCentroid(nodesj, xyz, jtype, elempt);
      Subtract(bpt, elempt, ne);
      CalcTriArea(&xyz[nodesi[0]*3], &xyz[nodesi[1]*3], &xyz[nodesi[2]*3], nb);
      Normalize(ne, ne);
      Normalize(nb, nb);
      dot = DotProduct(ne, nb);
      //reverse winding if normal doesn't face correct direction
      //that is, wind all face normals outward
      if(!(real(dot) > 0.0)){
	if(itype == TRI){
	  nodesitemp[0] = nodesi[2];	
	  nodesitemp[1] = nodesi[1];
	  nodesitemp[2] = nodesi[0];
	  elementi.Init(nodesitemp);
	} 
	else if(itype == QUAD){
	  nodesitemp[0] = nodesi[3];	
	  nodesitemp[1] = nodesi[2];
	  nodesitemp[2] = nodesi[1];
	  nodesitemp[3] = nodesi[0];
	  elementi.Init(nodesitemp);
	}
	else{
	  std::cerr << "WARNING!!! Unable to rewind element " << gelem << " Unknown type" << std::endl;
	  return (-1);
	}
	count++;
      }
      //std::cout << dot << std::endl;
    }
    gelem++;
  }
  std::cout << "MESH UTILITY: Rewinding " << count << " boundary elements" << std::endl;
  return (0);
}

//if reverse is true, reverse ordering is used
template <class Type>
Int Mesh<Type>::ReorderMeshCuthillMcKee(Int reverse)
{
  //a note about the STL choice here: deque is a double ended queue meaning
  //that access from both the end and beginning is constant time... list
  //was another alternative but it was deemed necessary to have random access
  //capability for list comparisons... vector is far too slow for large sets.
  using namespace std;
  Int i = 0;
  Int j, k;
  Int indx1, indx2;
  Int z;
  deque<Int> S;
  deque<IntInt> R;

  IntInt* Q = new IntInt[nnode];

  Int startnode = 0;
  S.clear();
  R.clear();

  cout << "MESH UTILITY: Reordering graph with Cuthill-McKee... " << std::endl;

  //init Q
  for(j = 0; j < nnode; j++){
    //here a is the index and b is the degree of freedom
    Q[j].a = j;
    Q[j].b = ipsp[j+1] - ipsp[j];
  }
  //choose first vertex in new ordering
  ordering[i] = startnode;
  //assign to seed
  z = ordering[i];
  Q[startnode].a = -1;

  while(!NegativeQueue(Q, nnode)){
    indx1 = ipsp[z];
    indx2 = ipsp[z+1];
    for(k = indx1; k < indx2; k++){
      //check if used
      if(psp[k] < nnode){
	if(Q[psp[k]].a >= 0){
	  R.push_back(Q[psp[k]]);
	  //set as used vertex
	  Q[psp[k]].a = -1;
	}
      }
    }
    sort(R.begin(), R.end(), DegreeCompare);
    while(!R.empty()){
      i++;
      ordering[i] = R.front().a;
      R.pop_front();
      S.push_back(ordering[i]);
    }
    if(!S.empty()){
      z = S.front();
      S.pop_front();
    }
    else{
      //Queue empty prematurely, pick a new seed
      //for now pick the first non-used seed in list Q
      for(j = 0; j < nnode; j++){
	if(Q[j].a >= 0){
	  z = j;
	  break;
	}
      }
    }
  }

  //reverse ordering
  if((bool)reverse){
    for(i = 0; i < nnode; i++){
      Q[i].a = ordering[i];
    }
    for(i = 0; i < nnode; i++){
      ordering[i] = Q[nnode-i-1].a;
    }
  }

#if 0
  //subsection of code prints bandwidth
  //point graphs out to file when enabled
  cout << "Writing bandwidth data..." << endl;
  ofstream fo1, fo2;
  Int kk;
  fo1.open("default_bandwidth.dat");
  fo2.open("CuthillMcKee_bandwidth.dat");
  for(i = 0; i < nnode; i++){
    //cout << ordering[i] << endl;
  }
  for(i = 0; i < nnode; i++){
    fo1 << i+1 << " " << i+1 << endl;
    indx1 = ipsp[i];
    indx2 = ipsp[i+1];
    for(k = indx1; k < indx2; k++){
      fo1 << i+1 << " " << psp[k]+1 << endl;
    }
  }
  for(i = 0; i < nnode; i++){
    j = ordering[i];
    fo2 << i+1 << " " << i+1 << endl;
    indx1 = ipsp[j];
    indx2 = ipsp[j+1];
    for(k = indx1; k < indx2; k++){
      fo2 << i+1 << " " ;
      for(kk = 0; kk < nnode; kk++){
	if(ordering[kk] == psp[k]){
	  fo2 << kk+1 << endl;
	  break;
	}
      }
    }
  }
  cout << "Finished writing bandwidth data" << endl;
  fo1.close();
  fo2.close();
#endif

  reordered = true;
 
  delete [] Q;
  return (0);
}

template <class Type>
Int Mesh<Type>::AllocateSolutionMemory()
{
  //create space for LSQ gradient coeff.
  s = new Type[(nnode+gnode)*6];
  sw = new Type[(nnode+gnode)*6];
  //allocate memory for mesh movement velocities
  nv = new Type[(nnode+gnode+nbnode)*3];
  xyzold = new Type[nnode*3];
  xyzoldm1 = new Type[nnode*3];

  //blank out the grid velocities since this only gets used with a moving grid
  MemBlank(nv, (nnode+gnode+nbnode)*3);

  solMemAlloc = 1;
  
  return(0);
}

template <class Type>
Int Mesh<Type>::IsInteriorNode(Int n) const
{
  return(n < nnode);
}

template <class Type>
Int Mesh<Type>::IsGhostNode(Int n) const
{
  return(n >= nnode && n < (nnode+gnode));
}

template <class Type>
Int Mesh<Type>::IsBoundaryNode(Int n) const
{
  return(n >= (nnode+gnode));
}

//destructor
template <class Type> 
Mesh<Type>::~Mesh()
{
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    delete *it;
  }
  
  delete [] nelem;
  delete [] mnode;
  delete [] extentsMax;
  delete [] extentsMin;
  if(meshInit == 1){
    delete [] xyz;
    delete [] xyz_base;
    delete [] vol;
    delete [] volold;
    delete [] vololdm1;
    delete [] cvstat;
    delete [] ordering;
    delete [] gNodeOwner;
    delete [] gNodeLocalId;
  }
  KillSolutionMemory();
  KillMapMemory();
  KillEdgesMemory();
}
 
template <class Type> 
Int Mesh<Type>::FindPointsWithFactag(Int** pts, Int factag)
{
  Int i, bedge, pt;
  Int indx, indx1, indx2;
  Int n = 0;
  Int memsize = 100;
  Int memjump = 50;
  *pts = new Int[memsize];

  for(i = 0; i < nsnodes; i++){
    pt = snodes[i];
    indx1 = ibesp[pt];
    indx2 = ibesp[pt+1];
    for(indx = indx1; indx < indx2; indx++){
      bedge = besp[indx];
      if(bedges[bedge].factag == factag){
	if(memsize <= n){
	  MemResize(pts, memsize, memsize+memjump);
	  memsize += memjump;
	}
	//make sure we don't accidentally add a ghost node to the list
	if(pt < nnode){
	  (*pts)[n] = pt;
	  n++;
	}
      }
    }
  }

  //do final memory resizing
  MemResize(pts, memsize, n);

  return n;
}

template <class Type> 
void Mesh<Type>::GetCoordsForPoints(Type* rxyz, Int* pts, Int n)
{
  Int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < 3; j++){
      rxyz[i*3 + j] = xyz[pts[i]*3 + j];
    }
  }
  return;
}

template <class Type>
Int Mesh<Type>::WriteCurrentCoords(std::string casename, Int timestep)
{
  Int err = 0;
  std::stringstream oss;
  oss.str("");
  oss << p->GetRank();
  Int np = p->GetNp();
  std::string h5filename = casename + "." + oss.str() + ".h5";
  hid_t h5out = HDF_OpenFile(h5filename, 1);
  if(h5out < 0){
    std::cerr << "Mesh::WriteParallelMesh() could not open file -- " << h5filename << std::endl;
    return (-1);
  }

  oss.str("");
  oss << timestep;
  std::string timestepstring = "timestep-" + oss.str() + "/";
  std::string directory = "/Mesh/" + timestepstring;

  //write local and ghost coords
  HDF_WriteArray(h5out, directory, "Nodal Coordinates", xyz, (nnode+gnode)*3);
  
  HDF_CloseFile(h5out);

  return err;
}

template <class Type>
Int Mesh<Type>::WriteParallelMesh(std::string casename)
{
  Int err = 0;
  std::stringstream oss;
  oss.str("");
  oss << p->GetRank();
  Int np = p->GetNp();
  std::string h5filename = casename + "." + oss.str() + ".h5";
  hid_t h5out = HDF_OpenFile(h5filename, 1);
  if(h5out < 0){
    std::cerr << "Mesh::WriteParallelMesh() could not open file -- " << h5filename << std::endl;
    return (-1);
  }

  //Write things which exist in the root directory
  HDF_WriteScalar(h5out, "/", "Number Of Processors", &np);
  Int reorder = reordered;
  HDF_WriteScalar(h5out, "/", "Reordered", &reorder);
  Int rescaled = scaled;
  HDF_WriteScalar(h5out, "/", "Rescaled", &rescaled);
  HDF_WriteScalar(h5out, "/", "Global Number Of Nodes", &gnnode);
  HDF_WriteScalar(h5out, "/", "Global Number Of Elements", &gnelem);
  HDF_WriteScalar(h5out, "/", "Number Of Factags", &nfactags);

  std::vector<Int> vint;
  std::vector<Int> vint2;
  //make a guess at the size needed to avoid repeat reallocations
  vint.reserve(lnelem*5);
  vint2.reserve(lnelem);
  
  std::string directory = "/Mesh/";

  HDF_WriteScalar(h5out, directory, "Number Of Triangles", &nelem[TRI]);
  HDF_WriteScalar(h5out, directory, "Number Of Quadrilaterals", &nelem[QUAD]);
  HDF_WriteScalar(h5out, directory, "Number Of Tetrahedron", &nelem[TET]);
  HDF_WriteScalar(h5out, directory, "Number Of Pyramids", &nelem[PYRAMID]);
  HDF_WriteScalar(h5out, directory, "Number Of Prisms", &nelem[PRISM]);
  HDF_WriteScalar(h5out, directory, "Number Of Hexahedron", &nelem[HEX]);

  //write number of local nodes
  HDF_WriteScalar(h5out, directory, "Number Of Local Nodes", &nnode);
  //write number of ghost nodes
  HDF_WriteScalar(h5out, directory, "Number Of Ghost Nodes", &gnode);
  //write ghost nodes owning process id
  HDF_WriteArray(h5out, directory, "Ghost Nodes Owning Process", gNodeOwner, gnode);
  //write local id of ghost nodes
  HDF_WriteArray(h5out, directory, "Ghost Nodes Local Id", gNodeLocalId, gnode);
  //write local and ghost coords
  HDF_WriteArray(h5out, directory, "Nodal Coordinates", xyz, (nnode+gnode)*3);
  //write elements
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int type = element.GetType();
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    vint.push_back(type);
    vint2.push_back(element.GetFactag());
    for(Int j = 0; j < nnodes; ++j){
      vint.push_back(nodes[j]);
    }
  }
  HDF_WriteArray(h5out, directory, "Element Data", &vint.front(), vint.size());
  HDF_WriteArray(h5out, directory, "Element Factags", &vint2.front(), vint2.size());

  HDF_CloseFile(h5out);

  return err;
}

template <class Type>
Int Mesh<Type>::ReadPartedMesh(std::string casename)
{
  Int i;

  Int npCheck;
  Int localElemCount;
  Int* factagTempLocal;
  std::ostringstream temposs;
  temposs.str("");
  temposs << p->GetRank();
  std::string h5filename = casename + "." + temposs.str() + ".h5";
  hid_t h5in = HDF_OpenFile(h5filename, 0);
  if(h5in < 0){
    std::cerr << "Mesh::ReadPartedMesh() could not open file -- " << h5filename << std::endl;
    return (-1);
  }

  std::string directoryBase =  "/Mesh/";
  std::cout << "PARTITION HDF I/O: Reading file --> " << h5filename << std::endl;
  
  //check to ensure number of processors is sane
  HDF_ReadScalar(h5in, "/", "Number Of Processors", &npCheck);
  if(npCheck != p->GetNp()){
    std::stringstream ss;
    ss << "WARNING: Number of processors does not match partition file!!\n";
    ss << "Running with " << p->GetNp() << " but file contains " << npCheck << "\n";
    HDF_CloseFile(h5in);
    Abort << ss.str();
    return (-1);
  }

  //read flag to tell if mesh has been previously reordered
  Int reorder;
  HDF_ReadScalar(h5in, "/", "Reordered", &reorder);
  reordered = reorder;

  //read flag to tell if mesh has been previously rescaled
  Int rescaled;
  HDF_ReadScalar(h5in, "/", "Rescaled", &rescaled);
  scaled = rescaled;

  //read total number of nodes in mesh
  HDF_ReadScalar(h5in, "/", "Global Number Of Nodes", &gnnode);
  //read total number of elems in mesh
  HDF_ReadScalar(h5in, "/", "Global Number Of Elements", &gnelem);
  //read number of factags in mesh
  HDF_ReadScalar(h5in, "/", "Number Of Factags", &nfactags);

  //read number of nodes in domain
  HDF_ReadScalar(h5in, directoryBase, "Number Of Local Nodes", &nnode);
  //read number of ghost nodes in domain
  HDF_ReadScalar(h5in, directoryBase, "Number Of Ghost Nodes", &gnode);
  //read number of each element type in domain
  HDF_ReadScalar(h5in, directoryBase, "Number Of Triangles", &nelem[TRI]);
  HDF_ReadScalar(h5in, directoryBase, "Number Of Quadrilaterals", &nelem[QUAD]);
  HDF_ReadScalar(h5in, directoryBase, "Number Of Tetrahedron", &nelem[TET]);
  HDF_ReadScalar(h5in, directoryBase, "Number Of Pyramids", &nelem[PYRAMID]);
  HDF_ReadScalar(h5in, directoryBase, "Number Of Prisms", &nelem[PRISM]);
  HDF_ReadScalar(h5in, directoryBase, "Number Of Hexahedron", &nelem[HEX]);

  //set values in mesh object
  localElemCount = 0;
  for(i = TRI; i <= HEX; i++){
    localElemCount += nelem[i];
  }
  bface = nelem[TRI] + nelem[QUAD];

  std::cout << "PARTITION HDF I/O: Number of local nodes " << nnode << std::endl;
  std::cout << "PARTITION HDF I/O: Number of ghost nodes " << gnode << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "PARTITION HDF I/O: Number of Hexes " << nelem[HEX] << std::endl;
  
  lnelem =  (nelem[TRI] + nelem[QUAD] + nelem[TET] + 
	     nelem[PYRAMID] + nelem[PRISM] + nelem[HEX]);
  
  MemInitMesh();

  //read in ghost nodes (owning domain)
  HDF_ReadArray(h5in, directoryBase, "Ghost Nodes Owning Process", &gNodeOwner, &gnode);
  //read in ghost nodes (local id on owning process)
  HDF_ReadArray(h5in, directoryBase, "Ghost Nodes Local Id", &gNodeLocalId, &gnode);

  //read in coords for local/ghost nodes
  //always read in Real values.. never complex
  Int ncoords = (nnode+gnode)*3;
  HDF_ReadArray(h5in, directoryBase, "Nodal Coordinates", &xyz, &ncoords);
  //we need the base xyz location for movement routines... copy over
  memcpy(xyz_base, xyz, sizeof(Type)*3*(nnode+gnode));

  factagTempLocal = new Int[localElemCount];
  //read in element factags
  HDF_ReadArray(h5in, directoryBase, "Element Factags", &factagTempLocal, &localElemCount);

  //read in element data
  Int elistSize = -1;
  Int* elementDataTemp = NULL;
  HDF_ReadArray(h5in, directoryBase, "Element Data", &elementDataTemp, &elistSize);
  elementDataTemp = new Int[elistSize];
  HDF_ReadArray(h5in, directoryBase, "Element Data", &elementDataTemp, &elistSize);

  //read in element data and put factags where they belong
  //they are stored as [type-e1, node1-e1, node2-e1, .... type-e2] in a linear array
  Int loc = 0;
  for(i = 0; i < localElemCount; i++){
    Int nodes[8]; 
    Int e = elementDataTemp[loc];
    loc++;
    for(Int j = 0; j < mnode[e]; j++){
      nodes[j] = elementDataTemp[loc];
      loc++;
    }
    Element<Type>* tempe = NULL;
    switch (e) {  
    case TRI :
      tempe = new Triangle<Type>;
      break;
    case QUAD :
      tempe = new Quadrilateral<Type>;
      break;
    case TET :
      tempe = new Tetrahedron<Type>;
      break;
    case PYRAMID :
      tempe = new Pyramid<Type>;
      break;
    case PRISM :
      tempe = new Prism<Type>;
      break;
    case HEX :
      tempe = new Hexahedron<Type>;
      break;
    default :
      std::cerr << "Type not defined in ReadPartedMesh()" << std::endl;
    }
    tempe->Init(nodes);
    tempe->SetFactag(factagTempLocal[i]);
    elementList.push_back(tempe);
  }
  std::cout << "PARTITION HDF I/O: File read successful!!" << std::endl;

  delete [] factagTempLocal;
  delete [] elementDataTemp;
  
  HDF_CloseFile(h5in);

  return (0);
}


template <class Type> template <class Type2>
Mesh<Type>::Mesh(const Mesh<Type2>& meshToCopy)
{
  //set all integer values to be equivalent
  lnelem = meshToCopy.GetNumLocalElem();
  nnode = meshToCopy.GetNumNodes();
  gnode = meshToCopy.GetNumParallelNodes();
  nbnode = meshToCopy.GetNumBoundaryNodes();
  gnelem = meshToCopy.GetNumGlobalElem();
  nedge = meshToCopy.GetNumEdges();
  ngedge = meshToCopy.GetNumParallelEdges();
  nbedge = meshToCopy.GetNumBoundaryEdges();
  bface = meshToCopy.bface;
  nfactags = meshToCopy.nfactags;
  nvnodes = meshToCopy.nvnodes;
  nsnodes = meshToCopy.nsnodes;

  DuplicateArray(&extentsMax, meshToCopy.extentsMax, 3);
  DuplicateArray(&extentsMin, meshToCopy.extentsMin, 3);

  //now deep copy all integer arrays
  DuplicateArray(&ordering, meshToCopy.ordering, nnode);
  for(Int i = 0; i < MAX_E_TYPES; ++i){
    mnode[i] = meshToCopy.GetNumElemNodes(i);
  }
  for(Int i = 0; i < MAX_E_TYPES; ++i){
    nelem[i] = meshToCopy.GetNumElem(i);
  }
  DuplicateArray(&gNodeOwner, meshToCopy.gNodeOwner, gnode);
  DuplicateArray(&gNodeLocalId, meshToCopy.gNodeLocalId, gnode);
  DuplicateArray(&ielsp, meshToCopy.ielsp, nnode+gnode+1);
  DuplicateArray(&elsp, meshToCopy.elsp, ielsp[nnode+gnode]);
  DuplicateArray(&iselsp, meshToCopy.iselsp, nnode+gnode+1);
  DuplicateArray(&selsp, meshToCopy.selsp, iselsp[nnode+gnode]);
  DuplicateArray(&ipsp, meshToCopy.ipsp, nnode+gnode+1);
  DuplicateArray(&psp, meshToCopy.psp, ipsp[nnode+gnode]);
  DuplicateArray(&iel2el, meshToCopy.iel2el, lnelem+1);
  DuplicateArray(&el2el, meshToCopy.el2el, iel2el[lnelem]);
  DuplicateArray(&isel2el, meshToCopy.isel2el, lnelem+1);
  Int size = 3*nelem[TRI] + 4*nelem[QUAD];
  DuplicateArray(&sel2el, meshToCopy.sel2el, size);
  DuplicateArray(&vnodes, meshToCopy.vnodes, nvnodes);
  DuplicateArray(&snodes, meshToCopy.snodes, nsnodes);
  DuplicateArray(&iesp, meshToCopy.iesp, nnode+gnode+1);
  DuplicateArray(&ibesp, meshToCopy.ibesp, nnode+gnode+1);
  DuplicateArray(&esp, meshToCopy.esp, iesp[nnode+gnode]);
  DuplicateArray(&besp, meshToCopy.besp, ibesp[nnode+gnode]);
  
  DuplicateArray(&xyz, meshToCopy.xyz, 3*nnode + 3*gnode);
  DuplicateArray(&xyzold, meshToCopy.xyzold, 3*nnode);
  DuplicateArray(&xyzoldm1, meshToCopy.xyzoldm1, 3*nnode);
  DuplicateArray(&xyz_base, meshToCopy.xyz_base, 3*nnode + 3*gnode);
  
  DuplicateArray(&vol, meshToCopy.vol, nnode);
  DuplicateArray(&volold, meshToCopy.volold, nnode);
  DuplicateArray(&vololdm1, meshToCopy.vololdm1, nnode);
  DuplicateArray(&cg, meshToCopy.cg, (nnode+gnode+nbnode)*3);

  DuplicateArray(&s, meshToCopy.s, (nnode+gnode)*6);
  DuplicateArray(&sw, meshToCopy.sw, (nnode+gnode)*6);
  DuplicateArray(&nv, meshToCopy.nv, (nnode+gnode+nbnode)*3);

  //allocate new edges
  edges = new Edges<RCmplx>[nedge];
  for(Int i = 0; i < nedge; i++){
    edges[i].n[0] = meshToCopy.edges[i].n[0];
    edges[i].n[1] = meshToCopy.edges[i].n[1];
    edges[i].wl = meshToCopy.edges[i].wl;
    edges[i].wr = meshToCopy.edges[i].wr;
    for(Int j = 0; j < 4; j++){
      edges[i].a[j] = meshToCopy.edges[i].a[j];
    }
  }
  bedges = new HalfEdges<RCmplx>[nbedge+ngedge];
  for(Int i = 0; i < (nbedge+ngedge); i++){
    bedges[i].n[0] = meshToCopy.bedges[i].n[0];
    bedges[i].n[1] = meshToCopy.bedges[i].n[1];
    bedges[i].elem = meshToCopy.bedges[i].elem;
    bedges[i].factag = meshToCopy.bedges[i].factag;
    for(Int j = 0; j < 4; j++){
      bedges[i].a[j] = meshToCopy.bedges[i].a[j];
    }
  }
 
  //allocate new control volume status flags
  //WARNING: these are not copied here but might should be
  cvstat = new CVStat[nnode];

  //copy element list
  for(typename std::vector<Element<Type2>*>::const_iterator it = meshToCopy.elementList.begin(); 
      it != meshToCopy.elementList.end(); ++it){
    //copy construct new elements to put onto the new mesh
    Element<Type2>& elemToCopy = **it;
    Int type = elemToCopy.GetType();
    Element<Type>* tmpElement;
    if(type == TRI){
      Triangle<Type2>& ecopy = static_cast<Triangle<Type2>&>(elemToCopy);
      tmpElement = new Triangle<Type>(ecopy);
    }
    else if(type == QUAD){
      Quadrilateral<Type2>& ecopy = static_cast<Quadrilateral<Type2>&>(elemToCopy);
      tmpElement = new Quadrilateral<Type>(ecopy);
    }
    else if(type == TET){
      Tetrahedron<Type2>& ecopy = static_cast<Tetrahedron<Type2>&>(elemToCopy);
      tmpElement = new Tetrahedron<Type>(ecopy);
    }
    else if(type == PYRAMID){
      Pyramid<Type2>& ecopy = static_cast<Pyramid<Type2>&>(elemToCopy);
     tmpElement = new Pyramid<Type>(ecopy);
    }
    else if(type == PRISM){
      Prism<Type2>& ecopy = static_cast<Prism<Type2>&>(elemToCopy);
      tmpElement = new Prism<Type>(ecopy);
    }
    else if(type == HEX){
      Hexahedron<Type2>& ecopy = static_cast<Hexahedron<Type2>&>(elemToCopy);
      tmpElement = new Hexahedron<Type>(ecopy);
    }
    else{
      Abort << "In Mesh(const Mesh&) and falling through element copy .. Buh Bye...";
    }
    
    elementList.push_back(tmpElement);
  }
  
  meshInit = true;
  mapsInit = true;
  mapsInitPsp = true;
  edgeInit = true;
  solMemAlloc = true;
  reordered = meshToCopy.IsReordered();
  scaled = meshToCopy.IsScaled();
}


template <class Type>
void Mesh<Type>::ScaleMesh(Type scaleFactor)
{
  //only do rescaling if it hasn't happened yet
  if(!scaled){
    std::cout << "MESH UTILITY: Scaling mesh by " << scaleFactor << std::endl;
    for(Int i = 0; i < 3*nnode + 3*gnode; i++){
      xyz[i] *= scaleFactor;
      xyz_base[i] *= scaleFactor;
    }
    scaled = true;
  }
}
