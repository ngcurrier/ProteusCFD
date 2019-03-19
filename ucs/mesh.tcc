#include "base64.h"
#include "math.h"
#include <algorithm>

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
Mesh<Type>::Mesh() :
  p(0), edges(0), bedges(0), cvstat(0),
  gNodeOwner(0), gNodeLocalId(0), cg(0), nv(0), s(0), sw(0),
  ordering(0), xyz(0), xyz_base(0), xyzold(0), xyzoldm1(0),
  vol(0), volold(0), vololdm1(0), psp(0), ipsp(0), besp(0),
  ibesp(0), nvnodes(0), vnodes(0), nsnodes(0), snodes(0),
  gnnode(0), nnode(0), gnode(0), nbnode(0), nelem(0), gnelem(0), lnelem(0),
  nedge(0), ngedge(0), nbedge(0), mnode(0), extentsMax(0), extentsMin(0),
  meshInit(false), mapsInit(false), mapsInitPsp(false), edgeInit(false),
  metricsCalc(false), solMemAlloc(false), scaled(false), reordered(false),
  esp(0), iesp(0), el2el(0), iel2el(0), sel2el(0), isel2el(0),
  elsp(0), ielsp(0), selsp(0), iselsp(0)
{

  nelem = new Int[MAX_E_TYPES];
  for(Int i = 0; i < MAX_E_TYPES; ++i){
    nelem[i] = 0;
  }
  
  //allocate space for extents
  extentsMax = new Type[3];
  extentsMin = new Type[3];

  //set number of vertices to expect for each elem type
  mnode = new Int[MAX_E_TYPES];
  mnode[TRI] = 3;
  mnode[QUAD] = 4;
  mnode[TET] = 4;
  mnode[PYRAMID] = 5;
  mnode[PRISM] = 6;
  mnode[HEX] = 8;
}

template <class Type>
Type const * Mesh<Type>::GetNodeCoords() const
{
  return xyz;
}

template <class Type>
Type* Mesh<Type>::GetNodeCoords()
{
  return xyz;
}

// This function returns the node coords prior to any movement routines being called
// and is useful for preventing roundoff accumulation during long movement simulations
template <class Type>
Type const * Mesh<Type>::GetNodeCoordsBase() const
{
  return xyz_base;
}

template <class Type>
void Mesh<Type>::SetParallelPointer(PObj<Type>* p_passed)
{
  p = p_passed;
}

template <class Type>
Type Mesh<Type>::GetVolumeTotal() const
{
  Type totVol = 0.0;
  for(Int i = 0; i < nnode; ++i){
    totVol += vol[i]; 
  }
  return totVol;
}

template <class Type>
Type const * Mesh<Type>::GetVolume() const
{
  return vol;
}

template <class Type>
Type const * Mesh<Type>::GetVolumeOld() const
{
  return volold;
}

template <class Type>
Type const * Mesh<Type>::GetVolumeOldM1() const
{
  return vololdm1;
}

template <class Type>
Int* Mesh<Type>::ElspBegin(Int ptid)
{
  Int indx1 = ielsp[ptid];
  return &(elsp[indx1]);
}

template <class Type>
Int* Mesh<Type>::ElspEnd(Int ptid)
{
  Int indx2 = ielsp[ptid+1];
  return &(elsp[indx2]);
}

template <class Type>
Int* Mesh<Type>::SelspBegin(Int ptid)
{
  Int indx1 = iselsp[ptid];
  return &(selsp[indx1]);
}

template <class Type>
Int* Mesh<Type>::SelspEnd(Int ptid)
{
  Int indx2 = iselsp[ptid+1];
  return &(selsp[indx2]);
}

template <class Type>
std::string Mesh<Type>::GetStringFromElemType(Int etype)
{
  switch (etype) {  
  case TRI :
    return "Triangle";
    break;
  case QUAD :
    return "Quadrilateral";
    break;
  case TET :
    return "Tetrahedron";
    break;
  case PYRAMID :
    return "Pyramid";
    break;
  case PRISM :
    return "Prism";
    break;
  case HEX :
    return "Hexahedron";
    break;
  default :
    std::cerr << "Type not defined in GetStringFromElemType()" << std::endl;
  }
  
}

template <class Type>
const Int* Mesh<Type>::ElspBegin(Int ptid) const
{
  Int indx1 = ielsp[ptid];
  return &(elsp[indx1]);
}

template <class Type>
const Int* Mesh<Type>::ElspEnd(Int ptid) const
{
  Int indx2 = ielsp[ptid+1];
  return &(elsp[indx2]);
}

template <class Type>
const Int* Mesh<Type>::SelspBegin(Int ptid) const
{
  Int indx1 = iselsp[ptid];
  return &(selsp[indx1]);
}

template <class Type>
const Int* Mesh<Type>::SelspEnd(Int ptid) const 
{
  Int indx2 = iselsp[ptid+1];
  return &(selsp[indx2]);
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
  xyzold = new Type[nnode*3 + 3*gnode];
  xyzoldm1 = new Type[nnode*3 + 3*gnode];
  vol = new Type[nnode];
  volold = new Type[nnode];
  vololdm1 = new Type[nnode];
  
  //initialize cv status flag objects
  cvstat = new CVStat[nnode];

  //this must be set here incase mesh reordering is not use
  //i.e. this must be a valid list as it is used in the buildMaps() function
  ordering = new Int[nnode];
  for(Int i = 0; i < nnode; ++i){
    ordering[i] = i;
  }

  //allocate parallel memory
  gNodeOwner = new Int[gnode];
  gNodeLocalId = new Int[gnode];

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
  //This function is valid to call after maps have been built
  //Therefore we have to attempt to clean up previous iterations
  //to avoid memory leaks
  KillSolutionMemory();
  KillMapMemory();
  KillEdgesMemory();
  
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

  //this has to be called last because we need to resolve the nbnode size
  //which is resolved from the BuildEdges() call
  AllocateSolutionMemory();
  
  mapsInit = true;
  mapsInitPsp = true;

  return(err);
}

//builds maps only up to psp() necessary for METIS interface
template <class Type>
Int Mesh<Type>::BuildMapsDecomp()
{
  //This function is valid to call after maps have been built
  //Therefore we have to attempt to clean up previous iterations
  //to avoid memory leaks
  KillSolutionMemory();
  KillMapMemory();
  KillEdgesMemory();

  int err = 0;
  err = BuildElsp();
  if(err != 0){
    return (err);
  }
  err = BuildPsp();
  if(err != 0){
    return (err);
  }
  
  mapsInitPsp = true;
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
  mapsInitPsp = false;
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
}

template <class Type>
void Mesh<Type>::KillSolutionMemory()
{
  if(solMemAlloc == 1){
    delete [] s;
    delete [] sw;
    delete [] nv;
  }
  solMemAlloc = false;
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
      if (node >= nnode){
	std::stringstream ss;
	ss << "MESH UTILITY: node found that is higher than node number in mesh\n";
	ss << "MESH UTILITY: node number is " << node << "\n";
	ss << "MESH UTILITY: number of nodes in mesh is " << nnode << "\n";
	ss << "MESH UTILITY: element containing is " << GetStringFromElemType(type);
	Abort << ss.str();
      }
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
      std::cerr << "MESH UTILITY: WARNING!!! Psp[] map failed sanity check on node " 
		<< i << std::endl;
      std::cerr << "MESH UTILITY: Node is likely an orphan, check mesh for validity and continuity in node numbering" << std::endl;
      
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
void Mesh<Type>::GetParallelTuple(Int node, Int& owningProcess, Int& localId) const
{
  if(!IsGhostNode(node)){
    owningProcess = -1;
    localId = -1;
    return;
  }
  else{
    node -= nnode;
    localId = gNodeLocalId[node];
    owningProcess = gNodeOwner[node];
    return;
  }
}

// computes the normal vector for a node based on the average of all neigboring
// surface element faces
// ptid - the node id for which to compute the normal vector
// normal - returned normal for the node based on neighbor faces
template <class Type>
void Mesh<Type>::GetNodeNeighborhoodNormal(Int ptid, std::vector<Type>& normal) const
{
  Int count = 0;
  std::vector<Type> avg(4,0.0);
  for(const Int* indx = SelspBegin(ptid); indx != SelspEnd(ptid); ++indx){
    Int selemid = *indx;
    std::vector<Type> inormal;
    elementList[selemid]->GetNormal(inormal, xyz);
    //all face normals are returned outward pointing, we want inward pointing. Flip them.
    avg[0] -= inormal[0]; //x
    avg[1] -= inormal[1]; //y
    avg[2] -= inormal[2]; //z
    avg[3] += inormal[3]; //actually face area
    count++;
  }
  avg[0] /= (Type)count;
  avg[1] /= (Type)count;
  avg[2] /= (Type)count;
  avg[3] /= (Type)count;
  normal = avg;
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
  std::cout << "MESH UTILITY: Calculating element areas and volumes" << std::endl;
  
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
  
  std::cout << "MESH UTILITY: Calculating face closure areas" << std::endl;
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
	  std::cout << "Negative volume detected in element " << jelem << " of type " << ETypeToName(jtype)
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
	  std::cout << "Negative volume detected in element " << jelem << " of type " << ETypeToName(jtype) 
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
  //initialize to zero for all nodes in 3 coordinate directions
  for(i = 0; i < nnode; i++){
    areasum[i*3 + 0] = 0.0;
    areasum[i*3 + 1] = 0.0;
    areasum[i*3 + 2] = 0.0;
  }
  //sum across all edges the 3 directional contributions
  //add from node 1 and subtract contributions from node 2 since areas are signed on the edge
  for(i = 0; i < nedge; i++){
    n1 = edges[i].n[0];
    n2 = edges[i].n[1];
    Subtract(&areasum[n1*3], edges[i].a, &areasum[n1*3]);
    Add(&areasum[n2*3], edges[i].a, &areasum[n2*3]);
  }
  //do the same thing with all boundary and ghost edges to ensure closure
  for(i = 0; i < ngedge+nbedge; i++){
    n1 = bedges[i].n[0];
    Subtract(&areasum[n1*3], bedges[i].a, &areasum[n1*3]);
  }
  //check that the dual of each node has faces which are closed, if not, show warning
  for(i = 0; i < nnode; i++){
    if(real(Magnitude(&areasum[i*3])) > real(smallnum)){
      fail = 1;
      std::cout << "WARNING!!! Dual CV based on node " << i << " not closed. Open by (" << areasum[i*3 + 0]
		<< ", " << areasum[i*3 + 1] << ", " << areasum[i*3 + 2] << ") -- position ("
		<< areasum[i*3 + 0] << " " << areasum[i*3 + 1] << " " 
		<< areasum[i*3 + 2] << ")" << std::endl;
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
		  << ETypeToName(itype) << " connected to " 
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
void Mesh<Type>::AllocateSolutionMemory()
{
  //create space for LSQ gradient coeff.
  s = new Type[(nnode+gnode)*6];
  sw = new Type[(nnode+gnode)*6];
  //allocate memory for mesh movement velocities
  nv = new Type[(nnode+gnode+nbnode)*3];

  //blank out the grid velocities since this only gets used with a moving grid
  MemBlank(nv, (nnode+gnode+nbnode)*3);

  solMemAlloc = 1;
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
    delete [] xyzold;
    delete [] xyzoldm1;
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

// returns (sorted) elements by id that are labeled with the appropriate factag
// elementIds - vector which will contain list of elements on return
// factag - factag of surface for which we'd like the elements
template <class Type> 
Int Mesh<Type>::FindSurfaceElementsWithFactag(std::vector<Int>& elementIds, Int factag)
{
  elementIds.resize(0);
  for(Int ielem = 0; ielem < elementList.size(); ++ielem){
    if(elementList[ielem]->GetFactag() == factag){
      elementIds.push_back(ielem);
    }
  }
  std::sort(elementIds.begin(), elementIds.end());
}


// returns sorted list of nodes by id that are labeled with the appropriate factag
// pts - unallocated pointer of Int which will contain list of nodes on return
// factag - factag of surface for which we are finding nodes
template <class Type> 
Int Mesh<Type>::FindPointsWithFactag(Int** pts, Int factag)
{
  Int i, bedge, pt;
  Int indx, indx1, indx2;
  Int n = 0;
  Int memsize = 100;
  Int memjump = 50;
  *pts = new Int[memsize];

  //make sure we don't add the node to the list twice
  Bool* used = new Bool[nnode];
  for(i = 0; i < nnode; ++i){
    used[i] = false;
  }
  
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
	if(pt < nnode && used[i] == false){
	  (*pts)[n] = pt;
	  n++;
	  used[i] = true;
	}
      }
    }
  }

  delete [] used;
  
  //do final memory resizing
  MemResize(pts, memsize, n);

  //sort list
  std::sort(*pts, *pts+n);

  return n;
}

//returns the largest numerical value of any factag
template <class Type>
Int Mesh<Type>::GetMaximumFactag(){
  Int max = -999;
  MPI_Datatype mpit;
  mpit = MPI_GetType(max);

  for(Int ielem = 0; ielem < elementList.size(); ++ielem){
    max = MAX(max, (elementList[ielem]->GetFactag()));
  }
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, mpit, MPI_MAX, MPI_COMM_WORLD);
  return max;
}

// returns coordinates for a list of points
// rxyz - array of coordinates, must be pre-allocated
// pts - array of points list by id
// n - number of points in pts list
template <class Type> 
void Mesh<Type>::GetCoordsForPoints(Type* rxyz, Int* pts, Int n)
{
  Int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < 3; j++){
      rxyz[i*3 + j] = xyz[pts[i]*3 + j];
    }
  }
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
  
  H5Fclose(h5out);

  return err;
}

template <class Type>
Int Mesh<Type>::GetNumElem() const
{
  return lnelem;
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

  H5Fclose(h5out);

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
    H5Fclose(h5in);
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
  
  H5Fclose(h5in);

  // do some basic sanity checks knowing that we need a field mesh to solve
  if(lnelem == 0){
    Abort << "WARNING: Mesh::ReadPartedMesh() found no local elements";
  }
  if(nelem[TET] +  nelem[PYRAMID] + nelem[PRISM] + nelem[HEX] == 0){
    Abort << "WARNING: Mesh::ReaPartedMesh() found no volume elements";
  }
  
  return (0);
}

// copy constructor for the mesh object
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

  extentsMax[0] = meshToCopy.GetMaxX();
  extentsMax[1] = meshToCopy.GetMaxY();
  extentsMax[2] = meshToCopy.GetMaxZ();
  extentsMin[0] = meshToCopy.GetMinX();
  extentsMin[1] = meshToCopy.GetMinY();
  extentsMin[2] = meshToCopy.GetMinZ();

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
  DuplicateArray(&xyzold, meshToCopy.xyzold, 3*nnode + 3*gnode);
  DuplicateArray(&xyzoldm1, meshToCopy.xyzoldm1, 3*nnode + 3*gnode);
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

// Append nodes to the current mesh, this assumes that the nodes adding are not ghost
// nodes but are on the mesh interior, does not fix elements, etc. and keep in mind that
// hanging nodes may wreak havoc on the core solver if not connected correctly
// see memInitMesh() for the things we have to adjust for
// numNewNodes - number of nodes we are appending
// xyz_new - list of new node coordinates  in flat array [x0,y0,z0, ..., zn,yn,zn)
template <class Type>
void Mesh<Type>::AppendNodes(Int numNewNodes, Type* xyz_new){
  Int numNodeOld = nnode;

  // the mesh needs nnode and gnode to be consistent with the size of the nnode array
  gnode = gnode;
  nnode = numNodeOld + numNewNodes;
  std::cout << "MESH UTILITY: appending " << numNewNodes << " new nodes to current " << numNodeOld << " nodes" << std::endl;
  
  // when we add nodes to the mesh we must adjust the following array sizes as well
  MemResize(&xyz, numNodeOld*3 + gnode*3, nnode*3 + gnode*3);
  MemResize(&xyz_base, numNodeOld*3 + gnode*3, nnode*3 + gnode*3);
  MemResize(&xyzold, numNodeOld*3 + gnode*3, nnode*3 + gnode*3);
  MemResize(&xyzoldm1, numNodeOld*3 + gnode*3, nnode*3 + gnode*3);

  MemResize(&vol, numNodeOld, nnode);
  MemResize(&volold, numNodeOld, nnode);
  MemResize(&vololdm1, numNodeOld, nnode);

  // don't use resize since not POD type
  delete [] cvstat;
  cvstat = new CVStat[nnode];
  
  MemResize(&ordering, numNodeOld, nnode);

  // we resized and copied across all the old data now we need to add in
  // the new data where appropriate

  Int j = 0;
  for(int i = numNodeOld*3; i < nnode*3; ++i){
    xyz[i] = xyz_new[j];
    xyz_base[i] = xyz[i];
    xyzold[i] = xyz[i];
    xyzoldm1[i] = xyz[i];
    ++j;
  }

  for(int i = 0; i < numNewNodes; ++i){
    vol[numNodeOld + i] = 0.0;
    volold[numNodeOld + i] = 0.0;
    vololdm1[numNodeOld + i] = 0.0;
  }

  // we just serialize the last nodes, they can be reordered later if we wish
  for(int i = 0; i < numNewNodes; ++i){
    ordering[numNodeOld + i] = numNodeOld + i;
  }
}


template <class Type>
Int Mesh<Type>::ReadSU2_Ascii(std::string filename)
{
  Element<Type>* tempe = NULL;
  std::ifstream fin;
  std::string line;
  Int ndim = 0;   //number of dimensions
  Int nnelem = 0;  //number of elements
  Int nbound = 0; //number of boundaries
    
  Int elemCounter = 0; //counts number of elements read
  Int nodeCount = 0; //counts number of nodes read
  Int boundElemCounter = 0; // counts boundary elements read
  Int boundCounter = 0; // counts boundary sections read
  Int elemType;
  Int nbelemTmp = 0; // temporary counter for each boundary section

  //
  // table converts from SU2 --> our format
  //
  Int translation[][8] = {
    {0,1,2}, // Triangle
    {0,1,2,3}, // Quad
    {0,1,2,3},  // Tet
    {0,1,2,3,4},  // Pyramid
    {0,1,2,3,4,5},  // Prism
    {0,1,2,3,4,5,6,7}  // Hex
  };
  
  //states for our state machine reader
  enum{
    stateReadNdim,
    stateReadNelem,
    stateReadElem,
    stateReadNpoints,
    stateReadPoints,
    stateReadNmark,
    stateReadMarkerTag,
    stateReadMarkerNelem,
    stateReadMarkerElem,
    stateExit,
    NSTATES
  };

  Int state = stateReadNdim;

  //the mesh is structured in the following order
  // 1) number of dimensions "NDIME=" -- we only accept 3d grids
  // 2) number of elements "NELEM="
  // 3) all elements (elemType, nodeIds, elementIndx)
  // 4) number of points "NPOIN="
  // 5) all point coordinates (x, y, z, id)
  // 6) number of defined boundaries "NMARK="
  // 7) defined boundaries (many sets defined as follows)
  //    * "MARKER_TAG= name"
  //    * number of elements on boundary "MARKER_ELEMS="
  //    * list of 2D boundary elements

  std::string ndimKey = "NDIME=";
  std::string nelemKey = "NELEM=";
  std::string npoinKey = "NPOIN=";
  std::string nmarkKey = "NMARK=";
  std::string markTagKey = "MARKER_TAG=";
  std::string markElemKey = "MARKER_ELEMS=";

  //NOTE: SU2 allows for comments with the % symbol
  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "Could not open file " << filename << " in SU2 ASCII I/O" << std::endl;
    return(-1);
  }

  std::cout << "SU2 ASCII I/O: Reading file --> " << filename << std::endl;

  //read each line one at a time, use a state machine to process data
  std::size_t loc;
  while(std::getline(fin, line)){
    if(line.size() == 0) continue; //skip blank lines
    if(line[0] == '%') continue; //skip comment lines
    std::stringstream ss;
    //use our state machine to do cool stuff with the datas...
    switch(state)
      {
      case stateReadNdim:
	loc = line.find(ndimKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NDIME=\" marker where expected\n";
	  return(1);
	}
	else{
	  ss.str(line.substr(loc+ndimKey.size()));
	  ss >> ndim;
	  std::cout << "SU2 ASCII I/O: Reading " << ndim << " dimensional mesh" << std::endl;
	  if(ndim != 3){
	    std::cerr << "SU2 ASCII I/O: File read only recognizes 3d meshes" << std::endl;
	    return(1);
	  }
	  else{
	    state = stateReadNelem;
	  }
	}
	break;
      case stateReadNelem:
	loc = line.find(nelemKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NELEM=\" marker where expected";
	  return(1);
	}
	else{
	  ss << line.substr(loc+nelemKey.size());
	  ss >> nnelem;
	  std::cout << "SU2 ASCII I/O: Reading " << nnelem << " volume elements" << std::endl;
	  state = stateReadElem;
	}
      	break;
      case stateReadElem:
	ss.str(line);
	ss >> elemType;
	if(elemType == SU2_LINE){
	  std::cerr << "SU2 ASCII I/O: Only 3d elements expected, found line";
	  return (-1);
	}
	else if(elemType == SU2_TRI){
	  std::cerr << "SU2 ASCII I/O: Only 3d elements expected, found triangle";
	  return (-1);
	}
	else if(elemType == SU2_QUAD){
	  std::cerr << "SU2 ASCII I/O: Only 3d elements expected, found quad";
	  return (-1);
	}
	else if(elemType == SU2_TET){
	  Int nodes[4];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[TET], TET, 0);
	  tempe = new Tetrahedron<Type>;
	  tempe->Init(nodes);
	  elementList.push_back(tempe);
	  nelem[TET]++;
	}
	else if(elemType == SU2_HEX){
	  Int nodes[8];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  ss >> nodes[4];
	  ss >> nodes[5];
	  ss >> nodes[6];
	  ss >> nodes[7];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
	  tempe = new Hexahedron<Type>;
	  tempe->Init(nodes);
	  elementList.push_back(tempe);
	  nelem[HEX]++;
	}
	else if(elemType == SU2_PRISM){
	  Int nodes[6];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  ss >> nodes[4];
	  ss >> nodes[5];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
	  tempe = new Prism<Type>;
	  tempe->Init(nodes);
	  elementList.push_back(tempe);
	  nelem[PRISM]++;
	}
	else if(elemType == SU2_PYRAMID){
	  Int nodes[5];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  ss >> nodes[4];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
	  tempe = new Pyramid<Type>;
	  tempe->Init(nodes);
	  elementList.push_back(tempe);
	  nelem[PYRAMID]++;
	}
	else{
	  std::cerr << "SU2 ASCII I/O: Cannot find element of type " << elemType << std::endl;
	  return(-1);
	}
	elemCounter++;
	if(elemCounter >= nnelem){
	  state = stateReadNpoints;
	}
	break;
      case stateReadNpoints:
	loc = line.find(npoinKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NPOIN=\" marker where expected\n";
	  return(1);
	}
	else{
	  ss.str(line.substr(loc+npoinKey.size()));
	  ss >> nnode;
	  std::cout << "SU2 ASCII I/O: Reading " << GetNumNodes() << " nodes" << std::endl;
	  state = stateReadPoints;
	  MemInitMesh();
	  std::cout << "SU2 ASCII I/O: Memory initialized" << std::endl;
	  std::cout << "SU2 ASCII I/O: Reading points " << std::endl;
	}
	break;
      case stateReadPoints:
	ss.str(line);
	{
	  Real xyzt[3];
	  ss >> xyzt[0];
	  ss >> xyzt[1];
	  ss >> xyzt[2];
	  Int id;
	  ss >> id;
	  xyz[nodeCount*3 + 0] = xyzt[0];
	  xyz[nodeCount*3 + 1] = xyzt[1];
	  xyz[nodeCount*3 + 2] = xyzt[2];
	  nodeCount++;
	}
	if(nodeCount >= nnode){
	  state = stateReadNmark;
	  std::cout << "SU2 ASCII I/O: Reading boundary markers" << std::endl;
	}
	break;
      case stateReadNmark:
	loc = line.find(nmarkKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NMARK=\" marker where expected" << std::endl;
	  return(1);
	}
	else{
	  ss.str(line.substr(loc+nmarkKey.size()));
	  ss >> nbound;
	  std::cout << "SU2 ASCII I/O: Reading " << nbound << " boundaries" << std::endl;
	  state = stateReadMarkerTag;
	}
	break;
      case stateReadMarkerTag:
	loc = line.find(markTagKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"MARKER_TAG=\" where expected\n";
	  return(1);
	}
	else{
	  std::string name;
	  ss.str(line.substr(loc+markTagKey.size()));
	  ss >> name;
	  std::cout << "SU2 ASCII I/O: Reading boundary " << name << std::endl;
	  state = stateReadMarkerNelem;
	}
	break;
      case stateReadMarkerNelem:
	loc = line.find(markElemKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"MARKER_ELEMS=\" where expected\n";
	  return(1);
	}
	else{
	  nbelemTmp = 0;
	  ss.str(line.substr(loc+markElemKey.size()));
	  ss >> nbelemTmp;
	  std::cout << "SU2 ASCII I/O: Reading " << nbelemTmp << " boundary elements " << std::endl;
	  state = stateReadMarkerElem;
	  boundElemCounter = 0;
	}
	break;
      case stateReadMarkerElem:
	ss.str(line);
	ss >> elemType;
	if(elemType == SU2_LINE){
	  std::cerr << "SU2 ASCII I/O: Only 3d elements expected, found line ";
	  return (-1);
	}
	else if(elemType == SU2_TRI){
	  Int nodes[3];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
	  tempe = new Triangle<Type>;
	  tempe->Init(nodes);
	  tempe->SetFactag(boundCounter+1);
	  elementList.push_back(tempe);
	  nelem[TRI]++;
	}
	else if(elemType == SU2_QUAD){
	  Int nodes[4];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  Int id;
	  ss >> id;
	  TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
	  tempe = new Quadrilateral<Type>;
	  tempe->Init(nodes);
	  tempe->SetFactag(boundCounter+1);
	  elementList.push_back(tempe);
	  nelem[QUAD]++;
	}
	else{
	  std::cerr << "SU2 ASCII I/O: Cannot find boundary element of type " << elemType << std::endl;
	  return(-1);
	}
	boundElemCounter++;
	if(boundElemCounter >= nbelemTmp){
	  state = stateReadMarkerTag;
	  boundCounter++;
	  if(boundCounter >= nbound){
	    state = stateExit;
	  }
	}
	break;
      default:
	break;
      }
    
  }

  nfactags = nbound;
  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }
  std::cout << "SU2 ASCII I/O: Number of nodes " << nnode << std::endl;
  std::cout << "SU2 ASCII I/O: Number of elements " << lnelem << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Hexes " << nelem[HEX] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of tagged boundaries " << nfactags << std::endl;

  fin.close();
  return(0);
}

// This function simply opens a gmsh file and then delegates the reading to the appropriate
// GMSH file version'd reader
template <class Type>
Int Mesh<Type>::ReadGMSH_Master(std::string filename)
{
  int ierr = 0;
  std::ifstream fin;
  std::string line;

  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "Could not open file " << filename << " in GMSH MASTER ASCII I/O" << std::endl;
    return(-1);
  }

  int majorVersion = 0;
  
  std::size_t loc;
  while(std::getline(fin, line)){
    if(line.size() == 0) continue; //skip blank lines
    if(line[0] == '%') continue; //skip comment lines
    //this is the mesh format indicator
    if(line.find("$MeshFormat") != std::string::npos){
      // read the next line which has the format in it
      // the first double value is the format id
      std::getline(fin, line);
      float version;
      std::stringstream ss(line);
      ss >> version;
      std::cout << "Found GMSH version: " << version << std::endl;
      majorVersion = floor(version);
    }
  }

  fin.close();
  if(majorVersion == 2 || majorVersion == 3){
    std::cout << "Opening file with GMSH format v2/3 style reader" << std::endl;
    ierr = ReadGMSH2_Ascii(filename);
  }
  else if(majorVersion == 4){
    std::cout << "Opening file with GMSH format v4 style reader" << std::endl;
    ierr = ReadGMSH4_Ascii(filename);
  }
  else{
    Abort << "Major version reader not available for this GMSH file type";
    
  }
  return ierr;
}

template <class Type>
Int Mesh<Type>::ReadGMSH2_Ascii(std::string filename)
{
  //Winding information (edgelist):
  // Triangle {0-1,1-2,2-0}
  // Quad {0-1,1-2,2-3,3-0}
  // Tet {0-3,3-1,1-0}, {3-1,1-2,2-3}, {1-0,0-2,2-1}, {0-3,3-2,2-0}
  // Prism {1-2,2-0,0-1}, {4-5,5-3,3-4}, {1-2,2-5,5-4,4-1}, {2-0,0-3,3-5,5-2}, {0-1,1-4,4-3,3-0}
  // Pyramid {3-0,0-1,1-2,2-3}, {1-2,2-4,4-1}, {0-1,1-4,4-0}, {2-3,3-4,4-2}, {3-0,0-4,4-3}
  // Hex {0-4,4-7,7-3,3-0}, {4-5,5-6,6-7,7-4}, {5-1,1-2,2-6,6-5}, {1-0,0-3,3-2,2-1}, {0-1,1-5,5-4,4-0}, {7-6,6-2,2-3,3-7}

  //
  // table converts from gmsh --> our format
  //
  Int translation[][8] = {
    {0,1,2}, // Triangle
    {0,1,2,3}, // Quad
    {0,1,2,3},  // Tet
    {0,1,2,3,4},  // Pyramid
    {0,1,2,3,4,5},  // Prism
    {0,1,2,3,4,5,6,7}  // Hex
  };

  Element<Type>* tempe = NULL;
  std::ifstream fin;
  std::string line;

  //states for our state machine reader
  enum s{
    None, /* this state is looking for new sections*/
    ReadMeshFormat,
    ReadNpoints,
    ReadPoints,
    ReadNelem,
    ReadElem,
    ReadNphysicalNames,
    ReadPhysicalNames,
    Exit,
    NSTATES
  };
  
  std::string versionTarget = "2.2";

  std::map<std::string, int> stateMap;
  stateMap["$MeshFormat"] = ReadMeshFormat;
  stateMap["$EndMeshFormat"] = None;
  stateMap["$ParametricNodes"] = ReadNpoints;
  stateMap["$EndParametricNodes"] = None;
  stateMap["$Nodes"] = ReadNpoints;
  stateMap["$EndNodes"] = None;
  stateMap["$Elements"] = ReadNelem;
  stateMap["$EndElements"] = None;
  stateMap["$PhysicalNames"] = ReadNphysicalNames;
  stateMap["$EndPhysicalNames"] = None;
  stateMap["$NodeData"] = None; //ignore this
  stateMap["$EndNodeData"] = None;
  stateMap["$ElementData"] = None; //ignore this
  stateMap["$EndElementData"] = None;
  
  Int state = None;
  Int nodesRead = 0; //counter to keep track of how many nodes have been read
  Int elemRead = 0;  //counter to keep track of how many elements have been read
  Int totalElems = 0; //tracks total elements expected
  
  std::cout << "GMSH ASCII I/O: Reading file --> " << filename << std::endl;

  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "Could not open file " << filename << " in GMSH ASCII I/O" << std::endl;
    return(-1);
  }

  //read each line one at a time, use a state machine to process data
  std::size_t loc;
  while(std::getline(fin, line)){
    line = trim(line);
    if(line.size() == 0) continue; //skip blank lines
    if(line[0] == '%') continue; //skip comment lines
    //this is a state change
    if(line[0] == '$'){
      std::map<std::string, int>::iterator it = stateMap.find(line);
      std::cout << "'" << line << "'" << std::endl;
      if(it != stateMap.end()){
	state = stateMap[line];
	continue;
      }
      else{
	std::cerr << "GMSH ASCII I/O: found bad state " << state << std::endl;
	std::cerr << "GMSH ASCII I/O: last line read " << line << std::endl;
	return(-1);
      }
    }
    std::stringstream ss;
    ss.str(line);
    //use our state machine to do cool stuff with the datas...
    Int numberOfTags = 0;
    std::vector<double> nodes;

    switch(state)
      {
      case None:
	//this just skips the line read and moves on
	break;
      case ReadMeshFormat:
	{
	  std::cout << "GMSH 2 ASCII I/O: Reading mesh format" << std::endl;
	  Real version;
	  Int tmp;
	  ss >> version;	
	  std::cout << "GMSH 2 ASCII I/O: Reading mesh format version " << version << std::endl;
	  ss >> tmp;
	  Int datasize;
	  ss >> datasize;
	  std::cout << "GMSH 2 ASCII I/O: Reading mesh datasize " << datasize << std::endl;
	}
	break;
      case ReadNpoints:
	ss >> nnode;
	std::cout << "GMSH 2 ASCII I/O: Number of points to read is " << nnode << std::endl;
	MemInitMesh();
	std::cout << "GMSH 2 ASCII I/O: Memory initialized" << std::endl;
	state = ReadPoints; //trip to next state for reading
	break;
      case ReadPoints:
	{
	  Int ptid;
	  ss >> ptid;
	  if(ptid-1 != nodesRead){
	    std::cerr << "GMSH 2 ASCII I/O: Warning points not ordered. Dying." << std::endl;
	    return -1;
	  }
	  ss >> xyz[nodesRead*3 + 0];
	  ss >> xyz[nodesRead*3 + 1];
	  ss >> xyz[nodesRead*3 + 2];
	  nodesRead++;
	}
	break;
      case ReadNelem:
	ss >> totalElems;
	std::cout << "GMSH 2 ASCII I/O: Number of elements to read is " << totalElems << std::endl;
	state = ReadElem; //trip to next state for reading
	break;
      case ReadElem:
	{
	  Int elemId = 0;
	  Int elemType = -1;
	  Int physicalId = 0;
	  Int elementaryId = 0;
	  // * first tag is the physical entity to which element belongs
	  // * second tag is the elementary geometric entity
	  // * third is the number of which partitions to which element belongs,
	  // * followed by the partition ids (negative implies ghost cells)
	  //     -- a zero tag is equivalent to no tag
	  // * format is elemId, elemType, number of tags, <tags> (physicalId, elementaryId, ... ) <list of nodes>
	  ss >> elemId;
	  ss >> elemType;
	  ss >> numberOfTags;
	  if(numberOfTags == 2){
	    ss >> physicalId;
	    ss >> elementaryId;
	  }
	  else{
	    std::cerr << "GMSH 2 ASCII I/O: number of tags inconsistent. Require two. One physical and one elementary tag id per element." << std::endl;
	    return (-1);
	  }
	  if(elemType == GMSH_LINE){
	    std::cerr << "GMSH 2 ASCII I/O: found line elements in file, these shouldn't be here. Ignoring.\n";
	  }
	  else if(elemType == GMSH_POINT){
	    std::cerr << "GMSH 2 ASCII I/O: found point elements in file, these shouldn't be here. Ignoring.\n";
	  }
	  else if(elemType == GMSH_TRI){
	    Int nodes[3];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
 	    tempe = new Triangle<Type>;
	    tempe->Init(nodes);
	    tempe->SetFactag(physicalId);
	    elementList.push_back(tempe);
	    nelem[TRI]++;
	  }
	  else if(elemType == GMSH_QUAD){
	    Int nodes[4];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
	    tempe = new Quadrilateral<Type>;
	    tempe->Init(nodes);
	    tempe->SetFactag(physicalId);
	    elementList.push_back(tempe);
	    nelem[QUAD]++;
	  }
	  else if(elemType == GMSH_TET){
	    Int nodes[4];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    TranslateWinding(nodes, translation, mnode[TET], TET, 0);
	    tempe = new Tetrahedron<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[TET]++;
	  }
	  else if(elemType == GMSH_HEX){
	    Int nodes[8];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    ss >> nodes[5];
	    ss >> nodes[6];
	    ss >> nodes[7];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    nodes[5]--;
	    nodes[6]--;
	    nodes[7]--;
	    TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
	    tempe = new Hexahedron<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[HEX]++;
	  }
	  else if(elemType == GMSH_PRISM){
	    Int nodes[6];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    ss >> nodes[5];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    nodes[5]--;
	    TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
	    tempe = new Prism<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[PRISM]++;
	  }
	  else if(elemType == GMSH_PYRAMID){
	    Int nodes[5];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
	    tempe = new Pyramid<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[PYRAMID]++;
	  }
	  else{
	    std::cerr << "GMSH 2 ASCII I/O: Element type " << elemType << " not valid" << std::endl;
	  }
	}
	break;
	
      default:
	break;
      }
  }
  fin.close();

  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }
  if(nelem[TRI] == 0 && nelem[QUAD] == 0){
    std::cerr << "GMSH 2 ASCII I/O: Warning no surface tri/quad elements found. Did you create physical groups on the surface for BCs?\n";
    return -1;
  }
  
  if(totalElems != lnelem){
    std::cerr << "GMSH 2 ASCII I/O: Number of elements expected does not match number read in" << std::endl;
    std::cerr << "GMSH 2 ASCII I/O: This sometimes happen if we ignore lower order elements like linear" << std::endl;
  }
  std::cout << "GMSH 2 ASCII I/O: Number of nodes " << nnode << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of elements " << lnelem << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of Hexes " << nelem[HEX] << std::endl;
  std::cout << "GMSH 2 ASCII I/O: Number of tagged boundaries " << nfactags << std::endl;
  
  return 0;
}

template <class Type>
Int Mesh<Type>::ReadGMSH4_Ascii(std::string filename)
{
  //Winding information (edgelist):
  // Triangle {0-1,1-2,2-0}
  // Quad {0-1,1-2,2-3,3-0}
  // Tet {0-3,3-1,1-0}, {3-1,1-2,2-3}, {1-0,0-2,2-1}, {0-3,3-2,2-0}
  // Prism {1-2,2-0,0-1}, {4-5,5-3,3-4}, {1-2,2-5,5-4,4-1}, {2-0,0-3,3-5,5-2}, {0-1,1-4,4-3,3-0}
  // Pyramid {3-0,0-1,1-2,2-3}, {1-2,2-4,4-1}, {0-1,1-4,4-0}, {2-3,3-4,4-2}, {3-0,0-4,4-3}
  // Hex {0-4,4-7,7-3,3-0}, {4-5,5-6,6-7,7-4}, {5-1,1-2,2-6,6-5}, {1-0,0-3,3-2,2-1}, {0-1,1-5,5-4,4-0}, {7-6,6-2,2-3,3-7}

  //
  // table converts from gmsh --> our format
  //
  Int translation[][8] = {
    {0,1,2}, // Triangle
    {0,1,2,3}, // Quad
    {0,1,2,3},  // Tet
    {0,1,2,3,4},  // Pyramid
    {0,1,2,3,4,5},  // Prism
    {0,1,2,3,4,5,6,7}  // Hex
  };

  Element<Type>* tempe = NULL;
  std::ifstream fin;
  std::string line;
  size_t lineCount = 0;
  
  //states for our state machine reader
  enum{
    None, /* this state is looking for new sections*/
    ReadMeshFormat,
    ReadNpoints,
    ReadPoints,
    ReadNelem,
    ReadElem,
    ReadNphysicalNames,
    ReadPhysicalNames,
    ReadEntityInfoNodes,
    ReadEntityInfoElems,
    ReadNEntities,
    ReadPointMapEntities,
    ReadCurveMapEntities,
    ReadSurfMapEntities,
    ReadVolMapEntities,
    Exit,
    NSTATES
  };
  
  std::map<std::string, int> stateMap;
  stateMap["$MeshFormat"] = ReadMeshFormat;
  stateMap["$EndMeshFormat"] = None;
  stateMap["$ParametricNodes"] = ReadNpoints;
  stateMap["$EndParametricNodes"] = None;
  stateMap["$Nodes"] = ReadNpoints;
  stateMap["$EndNodes"] = None;
  stateMap["$Elements"] = ReadNelem;
  stateMap["$EndElements"] = None;
  stateMap["$PhysicalNames"] = ReadNphysicalNames;
  stateMap["$EndPhysicalNames"] = None;
  stateMap["$NodeData"] = None; //ignore this
  stateMap["$EndNodeData"] = None;
  stateMap["$ElementData"] = None; //ignore this
  stateMap["$EndElementData"] = None;
  stateMap["$Entities"] = ReadNEntities;
  stateMap["$EndEntities"] = None;
  
  Int state = None;
  Int nodesRead = 0; //counter to keep track of how many nodes have been read
  Int elemRead = 0;  //counter to keep track of how many elements have been read
  Int totalElems = 0; //tracks total elements expected
  Int nEntityBlocks = 0; // tracks number of entities blocks
  Int nEntityNodes = 0; //local value for nodes in entity
  Int nEntityElems = 0; //local value for elems in entity
  Int entityNodesRead = 0; //counter for number of entity nodes read
  Int entityElemsRead = 0; //counter for number of entity elems read
  Int maxnodes = 0; //maximum node id that will be expected, used to reorder nodes for continuity

  
  std::cout << "GMSH 4 ASCII I/O: Reading file --> " << filename << std::endl;

  fin.open(filename.c_str());
  if(!fin.is_open()){
    std::cerr << "Could not open file " << filename << " in GMSH ASCII I/O" << std::endl;
    return(-1);
  }

  //read each line one at a time, use a state machine to process data
  std::size_t loc;
  std::vector<Bool> nodesReadBool;
  std::vector<Real> xyztmp;

  // GMSH files version > 3.0 do not have consecutive node numbering for reasons that can
  // only be described as really annoying.  We need to look for the $EndNodes tag and get
  // the real largest node id upfront to avoid segfaults and then rack and stack back
  // down the list to fill in the gaps
  int currPos = fin.tellg();
  int nowLine = fin.tellg();
  int prevLine = fin.tellg();
  std::string search = "$EndNodes";
  for(unsigned int curLine = 0; getline(fin, line); curLine++) {
    if (line.find(search) != std::string::npos) {
      std::string tmpLine;
      fin.seekg(prevLine);
      getline(fin, tmpLine);
      std::stringstream ss;
      std::vector<std::string> values = Tokenize(tmpLine, ' ');
      ss << values[0];
      ss >> maxnodes;
      break;
    }
    if(fin.eof()){
      Abort << "WARNING: file does not contain $EndNodes tag";
    }
    prevLine = nowLine;
    nowLine = fin.tellg();
  }
  // hop back to continue the reading
  fin.seekg(currPos);
  nodesReadBool.resize(maxnodes, false);
  xyztmp.resize(maxnodes*3, 0.0);
  line.clear();

  // First value in these maps is the gmsh identifier
  // Seocond value stored will be the physical identifier (i.e. boundary)
  Int nptentities, ncurveentities, nsurfentities, nvolentities;
  std::map<Int, Int> nodePhysIdMap;
  std::map<Int, Int> curvePhysIdMap;
  std::map<Int, Int> surfPhysIdMap;
  std::map<Int, Int> volPhysIdMap;
  while(std::getline(fin, line)){
    line = trim(line);
    lineCount++;
    if(line.size() == 0) continue; //skip blank lines
    if(line[0] == '%') continue; //skip comment lines
    //this is a state change
    if(line[0] == '$'){
      std::map<std::string, int>::iterator it = stateMap.find(line);
      if(it != stateMap.end()){
	state = stateMap[line];
	continue;
      }
      else{
	std::cerr << "GMSH 4 ASCII I/O: found bad state " << state << std::endl;
	std::cerr << "GMSH 4 ASCII I/O: last line read " << line << std::endl;
	return(-1);
      }
    }
    std::stringstream ss;
    ss.str(line);
    //use our state machine to do cool stuff with the datas...
    Int numberOfTags = 0;

    //entity information
    Int tagEntity;
    Int dimEntity;
    Int parametric;
    Int elemType;
    Int id;

    Real minx, miny, minz, maxx, maxy, maxz;
    Int nphysicalTags, ptag;
    
    switch(state)
      {
      case None:
	//this just skips the line read and moves on
	break;
      case ReadMeshFormat:
	{
	  std::cout << "GMSH 4 ASCII I/O: Reading mesh format" << std::endl;
	  Real version;
	  Int tmp;
	  ss >> version;	
	  std::cout << "GMSH 4 ASCII I/O: Reading mesh format version " << version << std::endl;
	  if (version != 4){
	    if( int(version) == 2){
	      return ReadGMSH2_Ascii(filename);
	    }
	    std::cerr << "GMSH 4 ASCII I/O: Expecting version 4.0. Quitting" << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: line # " << lineCount << std::endl;
	    return(-1);
	  }
	  ss >> tmp;
	  Int datasize;
	  ss >> datasize;
	  if(datasize != 8){
	    std::cerr << "GMSH 4 ASCII I/O: Data size expected is 8. Received " << datasize << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: line # " << lineCount << std::endl;
	  }
	}
	break;
      case ReadNEntities:
	{
	  std::cout << "GMSH 4 ASCII I/O: Reading entities header" << std::endl;
	  ss >> nptentities;
	  ss >> ncurveentities;
	  ss >> nsurfentities;
	  ss >> nvolentities;

	  std::cout << "GMSH 4 ASCII I/O: Number of point entities " << nptentities << std::endl;
	  std::cout << "GMSH 4 ASCII I/O: Number of curve entities " << ncurveentities << std::endl;
	  std::cout << "GMSH 4 ASCII I/O: Number of surface entities " << nsurfentities << std::endl;
	  std::cout << "GMSH 4 ASCII I/O: Number of volume entities " << nvolentities << std::endl;
	  state = ReadPointMapEntities;
	  break;
	}
      case ReadPointMapEntities:
	ss >> id;
	//points don't have a min/max coordinate, store in min only
	ss >> minx;
	ss >> miny;
	ss >> minz;
	ss >> nphysicalTags;
	ss >> ptag;
	nodePhysIdMap[id] = ptag;
	if(nodePhysIdMap.size() >= nptentities){
	  state = ReadCurveMapEntities;
	}
	break;
      case ReadCurveMapEntities:
	ss >> id;
	ss >> minx;
	ss >> miny;
	ss >> minz;
	ss >> maxx;
	ss >> maxy;
	ss >> maxz;
	ss >> nphysicalTags;
	ss >> ptag;
	curvePhysIdMap[id] = ptag;
	if(curvePhysIdMap.size() >= ncurveentities){
	  state = ReadSurfMapEntities;
	}
	break;
      case ReadSurfMapEntities:
	ss >> id;
	ss >> minx;
	ss >> miny;
	ss >> minz;
	ss >> maxx;
	ss >> maxy;
	ss >> maxz;
	ss >> nphysicalTags;
	ss >> ptag;
	surfPhysIdMap[id] = ptag;
	std::cout << id << "\tPTGA: " << ptag << std::endl;
	if(surfPhysIdMap.size() >= nsurfentities){
	  state = ReadVolMapEntities;
	}
	break;
      case ReadVolMapEntities:
	ss >> id;
	ss >> minx;
	ss >> miny;
	ss >> minz;
	ss >> maxx;
	ss >> maxy;
	ss >> maxz;
	ss >> nphysicalTags;
	ss >> ptag;
	volPhysIdMap[id] = ptag;
	if(volPhysIdMap.size() >= nvolentities){
	  state = None;
	}
	break;
      case ReadEntityInfoNodes:
	  ss >> tagEntity;
	  ss >> dimEntity;
	  ss >> parametric;
	  ss >> nEntityNodes;
	  std::cout << "GMSH 4 ASCII I/O: Reading entity " << tagEntity << " with " << nEntityNodes << " nodes" << std::endl;
 	  entityNodesRead = 0;
	  if(nEntityNodes != 0){
	    state = ReadPoints;
	  }
	  else{
	    state = ReadEntityInfoNodes;
	  }
	break;
      case ReadEntityInfoElems:
	  ss >> tagEntity;
	  ss >> dimEntity;
	  ss >> elemType;
	  ss >> nEntityElems;
	  Int physicalId;
	  if(surfPhysIdMap.count(tagEntity)){
	    physicalId = surfPhysIdMap[tagEntity];
	    //TODO: add support for volume tags to gmsh 4 reader
	  }
	  else{
	    physicalId = tagEntity;
	  }
	  std::cout << "GMSH 4 ASCII I/O: Reading entity " << tagEntity << " with " << nEntityElems << " elems - of type " << elemType << " with physical id "
		    << physicalId << std::endl;
	  entityElemsRead = 0;

	  if(nEntityElems != 0){
	    state = ReadElem;
	  }
	  else{
	    state = ReadEntityInfoElems;
	  }
	  break;
     case ReadNpoints:
       {
	ss >> nEntityBlocks;
	std::cout << "GMSH 4 ASCII I/O: Number of node entity blocks to read is " << nEntityBlocks << std::endl;
	ss >> nnode;
	std::cout << "GMSH 4 ASCII I/O: Number of points to read is " << nnode << std::endl;
	std::cout << "GMSH 4 ASCII I/O: Maximum node index is " << maxnodes << std::endl;
	if(nnode != maxnodes){
	  std::cout << "GMSH 4 ASCII I/O: Post read reordering will be required" << std::endl;
	}
	MemInitMesh();
	std::cout << "GMSH 4 ASCII I/O: Memory initialized" << std::endl;
	state = ReadEntityInfoNodes; //trip to next state for reading
	break;
       }
      case ReadPoints:
	{
	  Int ptid;
	  ss >> ptid;
	  if(nodesRead >= maxnodes || (ptid-1) >= maxnodes){
	    std::cerr << "GMSH 4 ASCII I/O: Total nodes read is " << nodesRead << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: Node read is " << ptid-1 << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: One of the above is greater than or equal to maximum " << nnode << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: line # " << lineCount << std::endl;
	    for(Int kk = 0; kk < nnode; ++kk){
	      if(nodesReadBool[kk] == false){
		std::cerr << "GMSH 4 ASCII I/O: Node missing " << kk << std::endl;
	      }
	    }
	    break;
	  }
	  nodesReadBool[ptid-1] = true;
	  ss >> xyztmp[(ptid-1)*3 + 0];
	  ss >> xyztmp[(ptid-1)*3 + 1];
	  ss >> xyztmp[(ptid-1)*3 + 2];
	  nodesRead++;
	  entityNodesRead++;
	  if(entityNodesRead >= nEntityNodes){
	    state = ReadEntityInfoNodes;
	  }
	}
	break;
      case ReadNelem:
	ss >> nEntityBlocks;
	std::cout << "GMSH 4 ASCII I/O: Number of element entity blocks to read is " << nEntityBlocks << std::endl;
	ss >> totalElems;
	std::cout << "GMSH 4 ASCII I/O: Number of elements to read is " << totalElems << std::endl;
	state = ReadEntityInfoElems; //trip to next state for reading
	break;
      case ReadElem:
	{
	  Int elemId = 0;
	  ss >> elemId;
	  Int physicalId;
	  if(surfPhysIdMap.count(tagEntity)){
	    physicalId = surfPhysIdMap[tagEntity];
	    //TODO: add support for volume tags to gmsh 4 reader
	  }
	  else{
	    physicalId = tagEntity;
	  }
	  if(elemType == GMSH_LINE){
	    std::cerr << "GMSH 4 ASCII I/O: found line elements in file, these shouldn't be here. Ignoring.\n";
	  }
	  else if(elemType == GMSH_POINT){
	    std::cerr << "GMSH 4 ASCII I/O: found point elements in file, these shouldn't be here. Ignoring.\n";
	  }
	  else if(elemType == GMSH_TRI){
	    Int nodes[3];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
	    tempe = new Triangle<Type>;
	    tempe->Init(nodes);
	    tempe->SetFactag(physicalId);
	    elementList.push_back(tempe);
	    nelem[TRI]++;
	  }
	  else if(elemType == GMSH_QUAD){
	    Int nodes[4];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
	    tempe = new Quadrilateral<Type>;
	    tempe->Init(nodes);
	    tempe->SetFactag(physicalId);
	    elementList.push_back(tempe);
	    nelem[QUAD]++;
	  }
	  else if(elemType == GMSH_TET){
	    Int nodes[4];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    TranslateWinding(nodes, translation, mnode[TET], TET, 0);
	    tempe = new Tetrahedron<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[TET]++;
	  }
	  else if(elemType == GMSH_HEX){
	    Int nodes[8];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    ss >> nodes[5];
	    ss >> nodes[6];
	    ss >> nodes[7];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    nodes[5]--;
	    nodes[6]--;
	    nodes[7]--;
	    TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
	    tempe = new Hexahedron<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[HEX]++;
	  }
	  else if(elemType == GMSH_PRISM){
	    Int nodes[6];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    ss >> nodes[5];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    nodes[5]--;
	    TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
	    tempe = new Prism<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[PRISM]++;
	  }
	  else if(elemType == GMSH_PYRAMID){
	    Int nodes[5];
	    ss >> nodes[0];
	    ss >> nodes[1];
	    ss >> nodes[2];
	    ss >> nodes[3];
	    ss >> nodes[4];
	    //gmsh nodes are 1 based, we are zero based
	    nodes[0]--;
	    nodes[1]--;
	    nodes[2]--;
	    nodes[3]--;
	    nodes[4]--;
	    TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
	    tempe = new Pyramid<Type>;
	    tempe->Init(nodes);
	    elementList.push_back(tempe);
	    nelem[PYRAMID]++;
	  }
	  else{
	    std::cerr << "GMSH 4 ASCII I/O: Element type " << elemType << " not valid" << std::endl;
	    std::cerr << "GMSH 4 ASCII I/O: line # " << lineCount << std::endl;
	  }
	  entityElemsRead++;
	  if(entityElemsRead >= nEntityElems){
	    state = ReadEntityInfoElems;
	  }
	}
	break;
	
      default:
	break;
      }
  }
  fin.close();

  // Do the element reordering if necessary
  std::cout << "GMSH 4 ASCII I/O: Performing node reordering step to ensure contiguous memory" << std::endl;
  std::vector<Int> missingNodes;
  Int countn = 0;
  for(Int kk = 0; kk < maxnodes; ++kk){
    if(nodesReadBool[kk] == false){
      missingNodes.push_back(kk);
    }
    else{
      xyz[countn*3 + 0] = xyztmp[kk*3 + 0];
      xyz[countn*3 + 1] = xyztmp[kk*3 + 1];
      xyz[countn*3 + 2] = xyztmp[kk*3 + 2];
      countn++;
    }
  }
  // compute the shift value for each node index
  std::vector<Int> nodeshift(maxnodes);
  Int loweridx = 0;
  Int shift = 0;
  for(Int i = 0; i < missingNodes.size(); ++i){
    Int skipidx = missingNodes[i];
    for(Int j = loweridx; j < skipidx; ++j){
      nodeshift[j] = shift;
    }
    shift++;
    loweridx = skipidx;
  }
  // get the number of nodes to shift by for every valid node that is present
  for(Int j = loweridx; j < maxnodes; ++j){
    nodeshift[j] = shift;
  }
  std::cout << "GMSH 4 ASCII I/O: Performing element rewinding to ensure correct node ids post reordering" << std::endl;
  Int maxelemnode = -1;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int etype = element.GetType();
    Int* nodes = NULL;
    // get node list and modify directly by precomputed shift
    Int ennodes = element.GetNodes(&nodes);
    for(Int j = 0; j < ennodes; ++j){
      Int n = nodes[j];
      nodes[j] = n - nodeshift[n];
      maxelemnode = MAX(maxelemnode, nodes[j]);
    }
  }  
  std::cout << "GMSH 4 ASCII I/O: Maximum node found in elements " << maxelemnode << std::endl;
  
  if(maxelemnode >= nnode){
    Abort << "GMSH 4 ASCII I/O: Elements found with node numbers that are not in range";
  }
  
  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }
  if(nelem[TRI] == 0 && nelem[QUAD] == 0){
    std::cerr << "GMSH 4 ASCII I/O: Warning no surface tri/quad elements found. Did you create physical groups on the surface for BCs?\n";
    return -1;
  }
  
  if(totalElems != lnelem){
    std::cerr << "GMSH 4 ASCII I/O: Number of elements expected does not match number read in" << std::endl;
    std::cerr << "GMSH 4 ASCII I/O: This sometimes happen if we ignore lower order elements like linear" << std::endl;
  }
  std::cout << "GMSH 4 ASCII I/O: Number of nodes " << nnode << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of elements " << lnelem << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of Hexes " << nelem[HEX] << std::endl;
  std::cout << "GMSH 4 ASCII I/O: Number of tagged boundaries " << nfactags << std::endl;

  return 0;
}

template <class Type>
void Mesh<Type>::UpdateElementCounts()
{
  nelem[TRI] = 0;
  nelem[QUAD] = 0;
  nelem[TET] = 0;
  nelem[PYRAMID] = 0;
  nelem[PRISM] = 0;
  nelem[HEX] = 0;
  nfactags = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int etype = element.GetType();
    nelem[etype]++;
    nfactags = MAX(nfactags, element.GetFactag());
  }

  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }

}

#ifdef _HAS_CGNS

template <class Type>
Int Mesh<Type>::ReadCGNS(std::string filename)
{
  int** isize = new int*[3];
  for(int i = 0; i < 3; ++i){
    isize[i] = new int[3];
  }
  Int error = GetCGNSSizes(filename, isize);

  MemInitMesh();

  CGNSreadCoordElem(filename, isize);

  for(int i = 0; i < 3; ++i){
    delete [] isize[i];
  }
  delete [] isize;
  
  return 0;
}

template <class Type>
Int Mesh<Type>::GetCGNSSizes(std::string filename, int** isize)
{
  int index_file;
  char zonename[100];

  cg_open(filename.c_str(), CG_MODE_READ, &index_file);

  int zone, nZones;
  int base = 1;
  if(cg_nzones(index_file, base, &nZones) != CG_OK)
    Abort << "readGridCGNS got nzone error";
  if(nZones != 1)
    Abort << "readGridCGNS got more than one zone";
  zone = 1;


  cg_zone_read(index_file, base, zone, zonename, isize[0]);

  nnode = isize[0][0];
  lnelem = isize[0][1];
  
  std::cout << "CGNS I/O: Number of nodes " << nnode << std::endl;
  std::cout << "CGNS I/O: Number of volume elements " << lnelem << std::endl;

  cg_close(index_file);

  return 0;
}

template <class Type>
void Mesh<Type>::CGNSreadCoordElem(std::string filename, int** isize)
{
  //This code borrows heavily from https://cgns.github.io/CGNS_docs_current/slides/VANDERWEIDE_tutorial.html
  
  int 			index_file,irmin,irmax,istart,iend,nbndry,iparent_flag,iparentdata;
  int 			*ielem, counter;
  ElementType_t 	itype;

  
  irmin = 1;
  irmax = isize[0][0];

  cg_open(filename.c_str(),CG_MODE_READ,&index_file);

  int nBases,base;
  if(cg_nbases(index_file, &nBases)!= CG_OK){
    Abort << "readGridCGNS got error";
    std::cerr << "CGNS error code " << cg_get_error() << std::endl;
  }
  if(nBases != 1){
    std::cerr << "CGNS error code " << cg_get_error() << std::endl;
    Abort << "readGridCGNS - reader assumes one base";
  }
  base = 1;

  char basename[100];
  int cellDim, physDim;
  if(cg_base_read(index_file, base, basename, &cellDim, &physDim) != CG_OK){
    std::cerr << "CGNS error code " << cg_get_error() << std::endl;
    Abort << "readGridCGNS got dimension error";
  }
  std::cout << "CGNS I/O: cell dimension " << cellDim << std::endl;
  std::cout << "CGNS I/O: physical dimension " << physDim << std::endl;
  std::cout << "CGNS I/O: reading base - " << std::string(basename) << std::endl;
  

  /* Read the number of zones in the grid. */
  int zone, nZones;
  if(cg_nzones(index_file, base, &nZones) != CG_OK)
    Abort << "readGridCGNS got nzone error";
  if(nZones != 1)
    Abort << "readGridCGNS got more than one zone";
  zone = 1;
  
  /* Check the zone type. This should be Unstructured. */
  ZoneType_t zoneType;
  if(cg_zone_type(index_file, base, zone, &zoneType) != CG_OK)
    Abort << "readGridCGNS got zone type error";
  if(zoneType != Unstructured)
    Abort << "readGridCGNS - Unstructured zone expected";
  else{
    std::cout << "CGNS I/O: reading zone - " << zone << std::endl;
  }


  //read coordinates
  Type* xtmp = new Type[isize[0][0]];
  cg_coord_read(index_file, base, zone, "CoordinateX", RealDouble, &irmin, &irmax, xtmp); //X
  for(int i = 0; i < nnode; ++i){
    xyz[i*3 + 0] = xtmp[i];
  }
  cg_coord_read(index_file, base, zone, "CoordinateY", RealDouble, &irmin, &irmax, xtmp); //Y
  for(int i = 0; i < nnode; ++i){
    xyz[i*3 + 1] = xtmp[i];
  }
  cg_coord_read(index_file, base, zone, "CoordinateZ", RealDouble, &irmin, &irmax, xtmp); //Z
  for(int i = 0; i < nnode; ++i){
    xyz[i*3 + 2] = xtmp[i];
  }
  delete [] xtmp;
  
  /* Determine the number of vertices and volume elements in this */
  /* We've already done this assuming one zone elsewhere*/

  /* Determine the number of sections for this zone. Note that */
  /* surface elements can be stored in a volume zone, but they */
  /* are NOT taken into account in the number obtained from */
  /* cg_zone_read. */

  int nSections;
  if(cg_nsections(index_file, base, zone, &nSections) != CG_OK){
    std::cerr << "CGNS error code " << cg_get_error() << std::endl;
    Abort << "readGridCGNS got error reading n-sections";
  }

  std::cout << "CGNS I/O: number of sections - " << nSections << std::endl;

  int ecount = 0;
  //read elements - note CGNS starts section numbering at 1 not zero
  for(int isec = 1; isec <= nSections; ++isec){
    char sectionname[100];
    if(cg_section_read(index_file, base, zone, isec, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag) != CG_OK){
      std::cerr << "CGNS error code " << cg_get_error() << std::endl;
      Abort << "readGridCGNS got error reading section";
    }
    std::cout << "CGNS I/O: reading section # " << isec << " - " << std::string(sectionname) << std::endl;
    cgsize_t elementDataSize;
    if(cg_ElementDataSize(index_file, base, zone, isec, &elementDataSize) != CG_OK){
      std::cerr << "CGNS error code " << cg_get_error() << std::endl;
      Abort << "readGridCGNS got error reading element data size";
    }
    std::cout << "CGNS I/O: element data size is - " << elementDataSize << std::endl;

    cgsize_t* conn = new cgsize_t[elementDataSize+1];

    if(cg_elements_read(index_file, base, zone, isec, conn, NULL) != CG_OK){
     std::cerr << "CGNS error code " << cg_get_error() << std::endl;
     Abort << "readGridCGNS got error reading connectivities";
    }

    int factagSection = isec; //using different terminology than Proteus but same effect
    switch (itype)
      {
      case TETRA_4:
	std::cout << "CGNS I/O: reading section of TETRA_4 elements" << std::endl;
	break;
      case PYRA_5:
	std::cout << "CGNS I/O: reading section of PYRA_5 elements" << std::endl;
	break;
      case PENTA_6:
	std::cout << "CGNS I/O: reading section of PENTA_6 elements" << std::endl;
	break;
      case HEXA_8:
	std::cout << "CGNS I/O: reading section of HEXA_8 elements" << std::endl;
	break;
      case TRI_3:
	std::cout << "CGNS I/O: reading section of TRI_3 elements" << std::endl;
	break;
      case QUAD_4:
	std::cout << "CGNS I/O: reading section of QUAD_4 elements" << std::endl;
	break;
      case MIXED:
	std::cout << "CGNS I/O: reading section of MIXED elements" << std::endl;
	break;
      default:
	// we should be aware that MIXED types are also support and contain
	// an extra integer to support identifying the element type
	Abort << "readGridCGNS Unsupported element encountered.";
	break;
      }

    //
    // table converts from cgns --> our format
    //
    Int translation[][8] = {
      {0,1,2}, // Triangle - same
      {0,1,2,3}, // Quad - same
      {0,1,2,3},  // Tet - same
      {0,1,2,3,4},  // Pyramid - same
      {0,1,2,3,4,5},  // Prism - same
      {0,1,2,3,4,5,6,7}  // Hex - same
    };
    
    int idx = 0;
    while(idx < elementDataSize){
      int type = itype;
      Int nodes[8];
      if(itype == MIXED){
	type = conn[idx]; //read in the first value which is element type
	idx++;
      }
      Element<Type>* tempe = NULL;
      switch (type){
      case TETRA_4:
	tempe = new Tetrahedron<Type>;
	for(int i = 0; i < 4; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[TET], TET, 0);
	idx+=4;
	nelem[TET]++;
	ecount++;
	break;
      case PYRA_5:
	tempe = new Pyramid<Type>;
	for(int i = 0; i < 5; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
	idx+=5;
	nelem[PYRAMID]++;
	ecount++;
	break;
      case PENTA_6:
	tempe = new Prism<Type>;
	for(int i = 0; i < 6; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
	idx+=6;
	nelem[PRISM]++;
	ecount++;
	break;
      case HEXA_8:
	tempe = new Hexahedron<Type>;
	for(int i = 0; i < 8; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
	idx+= 8;
	nelem[HEX]++;
	ecount++;
	break;
      case TRI_3:
	tempe = new Triangle<Type>;
	for(int i = 0; i < 3; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
	idx+=3;
	nelem[TRI]++;
	ecount++;
	break;
      case QUAD_4:
	tempe = new Quadrilateral<Type>;
	for(int i = 0; i < 4; ++i){
	  nodes[i] = conn[idx + i] - 1;
	}
	TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
	idx+=4;
	nelem[QUAD]++;
	ecount++;
	break;
      default:
	// we should be aware that MIXED types are also support and contain
	// an extra integer to support identifying the element type
	std::cerr << "Found element label type " << type << std::endl;
	std::cerr << "Type available are: " << TETRA_4 << " " << PYRA_5 << " " << PENTA_6 << " "
		  << HEXA_8 << " " << TRI_3 << " " << QUAD_4 << std::endl;
	Abort << "readGridCGNS Unsupported element encountered after parsing.";
	break;
      }
      tempe->Init(nodes);
      tempe->SetFactag(factagSection);
      elementList.push_back(tempe);
    }
    
    delete [] conn;
  }

  std::cout << "CGNS I/O: Read " << ecount << " elements" << std::endl;
  
  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }
  
  std::cout << "CGNS I/O: Number of nodes " << nnode << std::endl;
  std::cout << "CGNS I/O: Number of elements " << lnelem << std::endl;
  std::cout << "CGNS I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "CGNS I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "CGNS I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "CGNS I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "CGNS I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "CGNS I/O: Number of Hexes " << nelem[HEX] << std::endl;
  std::cout << "CGNS I/O: Number of tagged boundaries " << nSections << std::endl;

  
  cg_close(index_file);
}

template <class Type>
int Mesh<Type>::CGNSgetBCNum(char *name)
{
  int index_file,nbocos;

  cg_open(name,CG_MODE_READ,&index_file);
  cg_nbocos(index_file,1,1,&nbocos);
  cg_close(index_file);

  return nbocos-1;
}

template <class Type>
int Mesh<Type>::CGNSgetBCIndNum(char *name, int ib)
{
  int 			index_file,nbocos;
  int 			normalindex[3],ndataset;
  char 			boconame[100];
  BCType_t 		ibocotype;
  PointSetType_t	iptset;
  cgsize_t 		npts,normallistflag;
  DataType_t 		normaldatatype;

  cg_open(name,CG_MODE_READ,&index_file);
  cg_boco_info(index_file,1,1,ib,boconame,&ibocotype,
	       &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
  cg_close(index_file);

  return npts;
}

template <class Type>
void Mesh<Type>::CGNSreadBCConditions(char *name, int **bc)
{
  int 			index_file,nbocos,normallist;
  int 			normalindex[3],ndataset;
  char 			boconame[100];
  BCType_t 		ibocotype;
  PointSetType_t	iptset;
  cgsize_t 		npts,normallistflag;
  DataType_t 		normaldatatype;
  cgsize_t 		*ipnts;


  cg_open(name,CG_MODE_READ,&index_file);
  cg_nbocos(index_file,1,1,&nbocos);

  for (int ib = 2; ib <= nbocos; ib++)
    {
      cg_boco_info(index_file,1,1,ib,boconame,&ibocotype,
                   &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);

      ipnts    = new cgsize_t[npts];

      cg_boco_read(index_file,1,1,ib,ipnts,&normallist);
      for (int i=0; i < npts; i++)
	{
	  bc[ib-2][i] = ipnts[i];
	}

      delete ipnts;
    }

  cg_close(index_file);
}

#endif //end _HAS_CGNS

template <class Type>
Int Mesh<Type>::WriteCRUNCH_Ascii(std::string casename)
{
  Int i,j;
  Int etype;    
  Int tempnodes[8];

  std::ofstream fout;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);

  std::string filename = casename + ".crunch";
  
  //
  // table converts from crunch --> our format
  //
  Int translation[][8] = {
    {2,1,0}, // Triangle
    {3,2,1,0}, // Quad
    {0,1,3,2},  // Tet
    {0,1,2,3,4},  // Pyramid
    {0,3,4,1,5,2},  // Prism
    {0,1,3,2,4,5,7,6}  // Hex
  };


  fout.open(filename.c_str());
  if(!fout.is_open()){
    return(-1);
  }
  
  std::cout << "CRUNCH ASCII I/O: Writing file --> " << filename << std::endl;

  fout << "FIELDVIEW 2 5" << std::endl;  //write major_version minor_version

  //write header information
  fout << "CONSTANTS" << std::endl;
  fout << "0.0" << std::endl; //eltime
  fout << "1.0" << std::endl; //fsmach
  fout << "0.0" << std::endl; //alpha
  fout << "1.0" << std::endl; //re

  //write grids section
  fout << "GRIDS 1" << std::endl;

  //assume no surface IDs are present internally... so we name them
  fout << "Boundary Table " << nfactags << std::endl;
  for(i = 0; i < nfactags; i++){
    //use normal flag ("1") for all boundaries
    ss.clear();
    ss.str("");
    ss << "surface_";
    ss << i;
    fout << "1 " << ss.str() << std::endl;
  }
  
  fout << "VARIABLE NAMES 0" << std::endl;
  fout << "BOUNDARY VARIABLE NAMES 0" << std::endl;

  //fout << "BLOCK_ELEMENT_COUNT" << std::endl;
  
  fout << "NODES " << nnode << std::endl;

  //write nodes
  fout.setf(std::ios::scientific);
  fout.precision(15);
  for(i = 0; i < nnode; i++){
    fout << xyz[i*3 + 0] << " " << xyz[i*3 + 1] << " " << xyz[i*3 + 2] << std::endl;
  }
  
  //NOTE: CRUNCH nodes are 1 based so they need to be shifted
  //from zero base used internally

  //write boundary faces
  fout << "BOUNDARY FACES " << nelem[TRI]+nelem[QUAD] << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    etype = element.GetType();
    if(etype == TRI || etype == QUAD){
      Int* nodes = NULL;
      Int nnodes = element.GetNodes(&nodes);
      memcpy(tempnodes, nodes, sizeof(Int)*nnodes);
      TranslateWinding(tempnodes, translation, nnodes, etype, 1);
      //write tag and number of nodes for the boundary face
      fout << element.GetFactag() << " " << nnodes << " ";
      //write node ids for the face
      for(j = 0; j < nnodes; j++){
	fout << tempnodes[j]+1 << " ";
      }
      //endline
      fout << std::endl;
    }
  }
  
  //write volume elements
  fout << "ELEMENTS" << std::endl;
  Int crunch_type = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    etype = element.GetType();
    if(etype != TRI && etype != QUAD){
      switch(etype)
	{
	case TET: crunch_type = 1; break;
	case HEX: crunch_type = 2; break;
	case PRISM: crunch_type = 3; break;
	case PYRAMID: crunch_type = 4; break;	
	default:
	  std::cerr << "WARNING: non-identifiable volume element in WriteCRUNCH_Ascii()" 
		    << std::endl;
	  return (-1);
	}
      
      Int* nodes = NULL;
      Int nnodes = element.GetNodes(&nodes);
      
      //write element type and number of nodes
      fout << crunch_type << " " << nnodes << " ";
      //translate nodes to crunch winding
      memcpy(tempnodes, nodes, sizeof(Int)*nnodes);
      TranslateWinding(tempnodes, translation, nnodes, etype, 1);
      
      //write node ids
      for(j = 0; j < nnodes; j++){
	fout << tempnodes[j]+1 << " ";
      }
      //endline
      fout << std::endl;
    }
  }
  fout << "VARIABLES" << std::endl;

  fout.close();
  return (0);
}

template <class Type>
Int Mesh<Type>::WriteVTK_Ascii(std::string casename, std::vector<SolutionField<Type>*>& fields)
{
  Int j;
  std::ofstream fout;
  std::string filename = casename + ".vtk";

  std::cout << "VTK ASCII I/O: Writing solution file --> " << filename << std::endl;

  fout.open(filename.c_str());

  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << "Solver solution data" << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  fout << "POINTS " << nnode << " double" << std::endl;
  for(Int i = 0; i < nnode; i++){
    fout << xyz[i*3 + 0] << " " << xyz[i*3 + 1] << " " << xyz[i*3 + 2] << std::endl;
  }
  j = 0;
  for(Int e = TRI; e <= HEX; e++){
    j += nelem[e]*mnode[e];
  }
  j += lnelem;
  fout << "CELLS " << lnelem << " " << j << std::endl; 
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    fout << nnodes << " ";
    for(Int i = 0; i < nnodes; i++){
      fout << nodes[i] << " ";
    }
    fout << "\n";
  }

  fout << "CELL_TYPES " << lnelem << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int e = element.GetType();
    switch (e) {  
    case TRI :
      fout << "5" << std::endl;
      break;
    case QUAD :
      fout << "9" << std::endl;
      break;
    case TET :
      fout << "10" << std::endl;
      break;
    case PYRAMID :
      fout << "14" << std::endl;
      break;
    case PRISM :
      fout << "13" << std::endl;
      break;
    case HEX :
      fout << "12" << std::endl;
      break;
    default :
      std::cerr << "Type not defined WriteVTK_Ascii() -- type " << e << std::endl;
    }
  }
  
  fout << "CELL_DATA " << lnelem << std::endl;
  fout << "SCALARS BoundaryIds int " << 1 << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int e = element.GetFactag();
    fout << e << "\n" ;
  }
  if(fields.size() != 0){
    fout << "POINT_DATA " << nnode << std::endl;
  }
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type>& field = **it;
    Int neqn = field.GetNdof();
    Real* q = field.GetData(FIELDS::STATE_NONE);
    for(j = 0; j < neqn; j++){
      if(field.DofIsVector(j)){
	fout << "VECTORS " << field.GetDofName(j) << " double " << std::endl;
	for(Int i = 0; i < nnode; i++){
	  fout << q[i*neqn + j] << " ";
	  fout << q[i*neqn + j+1] << " ";
	  fout << q[i*neqn + j+2] << std::endl;;
	}
	j+=2;
      }
      else{
	fout << "SCALARS " << field.GetDofName(j) << " double " << 1 << std::endl;
	fout << "LOOKUP_TABLE default" << std::endl;
	for(Int i = 0; i < nnode; i++){
	  fout << q[i*neqn + j] << "\n";
	}
      }
    }
  }
  fout << std::endl;

  fout.close();

  std::cout << "VTK ASCII I/O: File write successful!!" << std::endl;


  return (0);
}


template <class Type>
Int Mesh<Type>::WriteGridXDMF(PObj<Type> &p, std::string filebase, std::string meshname)
{
  Int i;
  Int err = 0;
  Int rank = p.GetRank();
  Int np = p.GetNp();
  Int nelem[np];
  Int nnode[np];
  Int elemSize[np];
  std::stringstream ss;
  std::string datapath;
  std::string h5filename = filebase + ".h5";
  std::string filename;
  
  //element types in XDMF are in the following order - 1 indexed
  // Polyvertex - a group of unconnected points
  // Polyline - a group of line segments
  // Polygon
  // Triangle
  // Quadrilateral
  // Tetrahedron
  // Pyramid
  // Wedge
  // Hexahedron 
  
  nnode[0] = nnode+gnode;

  //testing
  nelem[0] = nnode[0];
  elemSize[0] = nelem[0]*2;
  Int* elements = new Int[nnode[0]*2];
  for(i = 0; i < nnode[0]; i++){
    elements[i*2 + 0] = 1;
    elements[i*2 + 1] = i;
  }

  //parallel sync number of elements of each type from each process

  //only allow root to write xml file
  if(rank == 0){
    TiXmlDocument doc;
    
    TiXmlElement* xdmf = new TiXmlElement("Xdmf");
    xdmf->SetAttribute("Version", 2.0);
    doc.LinkEndChild(xdmf);
    
    TiXmlElement* domain = new TiXmlElement("Domain");
    xdmf->LinkEndChild(domain);

    for(i = 0; i < np; i++){
      ss.str("");
      ss << i;
      std::string meshNumbered = meshname + ss.str();
      TiXmlElement* grid = new TiXmlElement("Grid");
      grid->SetAttribute("Name", meshNumbered.c_str());
      grid->SetAttribute("GridType", "Uniform");
      domain->LinkEndChild(grid);
      
      TiXmlElement* topology = new TiXmlElement("Topology");
      topology->SetAttribute("TopologyType", "Mixed");
      topology->SetAttribute("NumberOfElements", nelem[i]);
      grid->LinkEndChild(topology);
      TiXmlElement* elementsData = new TiXmlElement("DataStructure");
      elementsData->SetAttribute("Dimensions", elemSize[i]);
      elementsData->SetAttribute("NumberType", "Int");
      elementsData->SetAttribute("Format", "HDF");
      datapath = "/DOM-" + ss.str() + "/mesh/connectivity";
      datapath = h5filename + ":" + datapath;
      elementsData->LinkEndChild( new TiXmlText(datapath.c_str()));
      topology->LinkEndChild(elementsData);
			        
      TiXmlElement* geometry = new TiXmlElement("Geometry");
      geometry->SetAttribute("GeometryType", "XYZ");
      grid->LinkEndChild(geometry);
      TiXmlElement* coordsData = new TiXmlElement("DataStructure");
      coordsData->SetAttribute("Dimensions", 3*nnode[i]);
      coordsData->SetAttribute("NumberType", "Float");
      coordsData->SetAttribute("Precision", 8);
      coordsData->SetAttribute("Format", "HDF");
      datapath = "/DOM-" + ss.str() + "/mesh/coords";
      datapath = h5filename + ":" + datapath;
      coordsData->LinkEndChild( new TiXmlText(datapath.c_str()));
      geometry->LinkEndChild(coordsData);
    }

    //dump_to_stdout(&doc);
    filename = filebase + ".xmf";
    doc.SaveFile(filename.c_str());
  }

  //Now do the parallel writing part to the HDF file
  hid_t fileId;
  fileId = HDF_OpenFile(h5filename, 1);

  std::string directory;
  std::string dataname;
  
  //write connectivity
  ss.str("");
  ss << rank;
  directory = "/DOM-" + ss.str() + "/mesh/";
  dataname = "connectivity";
  std::cout << "Writing to : " << directory+dataname << std::endl;
  HDF_WriteArray(fileId, directory, dataname, elements, nnode[0]*2);
    

  //write xyz coords
  ss.str("");
  ss << rank;
  directory = "/DOM-" + ss.str() + "/mesh/";
  dataname = "coords";
  std::cout << "Writing to : " << directory+dataname << std::endl;
  HDF_WriteArray(fileId, directory, dataname, xyz, (nnode+gnode)*3);
						     
  H5Fclose(fileId);

  delete [] elements;

  return err;

}

template <class Type>
Int Mesh<Type>::ReadUGRID_Ascii(std::string filename)
{
  std::ifstream fin;
  Int i, j;
  Bool tag_shift = 0;
  Int nodes[8];
  Element<Type>* tempe = NULL;

  //
  // table converts from ugrid --> our format
  //
  Int translation[][8] = {
    {2,1,0}, // Triangle
    {3,2,1,0}, // Quad
    {0,1,2,3},  // Tet
    {0,3,4,1,2},  // Pyramid
    {0,1,2,3,4,5},  // Prism
    {0,1,2,3,4,5,6,7}  // Hex
  };

  fin.open(filename.c_str());
  if(!fin.is_open()){
    return(-1);
  }
  
  std::cout << "UGRID ASCII I/O: Reading file --> " << filename << std::endl;

  fin >> nnode >> nelem[TRI] >> nelem[QUAD] >> nelem[TET];
  fin >> nelem[PYRAMID] >> nelem[PRISM] >> nelem[HEX];

  std::cout << "UGRID ASCII I/O: Number of nodes " << nnode << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Hexes " << nelem[HEX] << std::endl;

  lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    lnelem += nelem[e];
  }

  //call init routine for mesh
  MemInitMesh();
  std::cout << "UGRID ASCII I/O: Memory initialized" << std::endl;

  std::cout << "UGRID ASCII I/O: Reading coordinates" << std::endl;

  //read in nodes
  for(i = 0; i < nnode; i++){
    for(j = 0; j < 3; j++){
      fin >> xyz[3*i + j];
    } 
  }
  std::cout << "UGRID ASCII I/O: Reading boundary faces" << std::endl;


  //NOTE: UGRID nodes are 1 based so they need to be shifted
  //to be zero based... ce la vi
  
  //read in Tris
  for(i = 0; i < nelem[TRI]; i++){
    for(j = 0; j < 3; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
    tempe = new Triangle<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }
  //read in Quads
  //UGRID winds these backwards
  for(i = 0; i < nelem[QUAD]; i++){
    for(j = 0; j < 4; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
    tempe = new Quadrilateral<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }
  
  //read in boundary tags for tris
  nfactags = 0;
  Int ielem = 0;
  for(i = 0; i < nelem[TRI]; i++){
    Int temp;
    fin >> temp;
    if(temp < 0){
      std::cout << "WARNING!! negative tags found.. DO NOT USE" << std::endl;
    }
    if(temp == 0){
      tag_shift = 1;
    }
    elementList[ielem]->SetFactag(abs(temp));
    nfactags = MAX(nfactags, elementList[ielem]->GetFactag());
    ielem++;
  }
  //read in boundary tags for quads
  for(i = 0; i < nelem[QUAD]; i++){
    Int temp;
    fin >> temp;
    if(temp < 0){
      std::cout << "WARNING!! negative tags found.. DO NOT USE" << std::endl;
    }
    if(temp == 0){
      tag_shift = 1;
    }
    elementList[ielem]->SetFactag(abs(temp));
    nfactags = MAX(nfactags, elementList[ielem]->GetFactag());
    ielem++;
  }
  

  //flip all tags for boundaries to negative
  //and if a zero tag was used shift all tags
  //away -- zero is used for (no tag) condition
  if(tag_shift == 1){
    std::cout << "\n\n";
    std::cout << "UGRID ASCII I/O: WARNING!!! ZERO TAGS RESERVED... DO NOT USE" << std::endl;
    std::cout << "UGRID ASCII I/O: WARNING!!! ATTEMPTING TO SHIFT TAG IDS BY ONE" << std::endl;
    std::cout << "UGRID ASCII I/O: WARNING!!! TAGS WILL BE WRITTEN OUT SHIFTED" << std::endl;
    std::cout << "\n\n";
    for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
      Element<Type>& element = **it;
      Int etype = element.GetType();
      if(etype == TRI || etype == QUAD){
	Int factag = element.GetFactag();
	element.SetFactag(factag + 1);
      }
    }
  }
  std::cout << "UGRID ASCII I/O: Reading volume mesh" << std::endl;

  //read in Tets
  for(i = 0; i < nelem[TET]; i++){
    for(j = 0; j < 4; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[TET], TET, 0);
    tempe = new Tetrahedron<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }
  //read in Pyramids
  //UGRID winds these weird... we don't use their winding
  for(i = 0; i < nelem[PYRAMID]; i++){
    for(j = 0; j < 5; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
    tempe = new Pyramid<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }
  //read in Prisms
  for(i = 0; i < nelem[PRISM]; i++){
    for(j = 0; j < 6; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
    tempe = new Prism<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }
  //read in Hexes
  for(i = 0; i < nelem[HEX]; i++){
    for(j = 0; j < 8; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
    tempe = new Hexahedron<Type>;
    tempe->Init(nodes);
    elementList.push_back(tempe);
  }

  std::cout << "UGRID ASCII I/O: File read successful!!" << std::endl;

  fin.close();
  return(0);
}

template <class Type>
Int Mesh<Type>::ReadCRUNCH_Ascii(std::string filename)
{
  Element<Type>* tempe = NULL;
  std::ifstream fin;
  Int i, j;
  Int count;
  Bool tag_shift = 0;
  Int indx1, indx2, indx3, indx4;

  std::string trash;
  std::string elemline = "elements";
  int elemloc = 0;
  std::string btableline = "boundary table";
  std::string bline = "boundary";
  std::string bfaceline = "boundary faces";
  Int bfaceloc = 0;
  std::string nodesline = "nodes";
  Int nodesloc = 0;
  
  Int numbfaces = 0;
  Int numvolelem = 0;

  Int line;
  size_t loc;

  Int tag;
  Int temp;
  Int npts;
  Int elemtype;
  Int nodes[8];

  enum
  {
    Cblank,
    CTet,
    CHex,
    CPrism,
    CPyramid
  };

  //
  // table converts from crunch --> our format
  //
  Int translation[][8] = {
    {2,1,0}, // Triangle
    {3,2,1,0}, // Quad
    {0,1,3,2},  // Tet
    {0,1,2,3,4},  // Pyramid
    {0,3,4,1,5,2},  // Prism
    {0,1,3,2,4,5,7,6}  // Hex
  };

  fin.open(filename.c_str());
  if(!fin.is_open()){
    return(-1);
  }
  
  std::cout << "CRUNCH ASCII I/O: Reading file --> " << filename << std::endl;

  //skip through everything in the header looking for "Boundary Table"
  //we just don't care about the rest at this point, don't look past 50 items
  count = 0;
  do{
    trash = "";
    line = fin.tellg();
    fin >> trash;
    trash = UpperToLower(trash);
    loc = trash.find(bline);
    count++;
  }while((loc == std::string::npos) && (count < 50));
  
  if(loc == std::string::npos){
    std::cerr << "CRUNCH ASCII I/O: BOUNDARY TABLE HEADER NOT IN EXPECTED LOCATION -- "
	      << "remove comments from file (i.e. lines with !)" << std::endl;
    return(1); 
  }

  loc += (btableline.size() + line + 1);
  fin.seekg(loc);
  fin >> nfactags;
  getline(fin, trash);
  //skip boundary table lines
  for(i = 0; i < nfactags; i++){
    getline(fin, trash);
  }
  trash.clear();

  //start looking for "Nodes"
  count = 0;
  do{
    trash = "";
    line = fin.tellg();
    fin >> trash;
    trash = UpperToLower(trash);
    loc = trash.find(nodesline);
    count++;
  }while((loc == std::string::npos) && (count < 50));

  if(loc == std::string::npos){
    std::cerr << "CRUNCH ASCII I/O: NODES HEADER NOT IN EXPECTED LOCATION -- "
	      << "remove comments from file (i.e. lines with !)" << std::endl;
    return(1);
  }

  loc += (nodesline.size() + line + 1);
  fin.seekg(loc);
  fin >> nnode;
  getline(fin, trash);
  nodesloc = fin.tellg();
  //burn through node lines... we'll read them later
  for(i = 0; i < nnode; i++){
    getline(fin, trash);
  }
  trash.clear();
  line = fin.tellg();
  getline(fin, trash);
  //convert to lowercase
  trash = UpperToLower(trash);
  loc = trash.find(bfaceline);
  if(loc == std::string::npos){
    std::cerr << "CRUNCH ASCII I/O: BOUNDARY FACES HEADER NOT IN EXPECTED LOCATION -- "
	      << "remove comments from file (i.e. lines with !)" << std::endl;
    return(1);
  }
  loc += (bfaceline.size() + line);
  fin.seekg(loc);
  fin >> numbfaces;
  getline(fin, trash);
  bfaceloc = fin.tellg();
  nelem[TRI] = 0;
  nelem[QUAD] = 0;
  for(i = 0; i < numbfaces; i++){
    fin >> temp >> npts;
    getline(fin, trash);
    if(npts == 3){
      nelem[TRI]++;
    }
    else if(npts == 4){
      nelem[QUAD]++;
    }
    else{
      std::cerr << "CRUNCH ASCII I/O: Unidentified number of points in boundary face " 
		<< npts << std::endl;
      return(1);
    }
  }
  trash.clear();
  line = fin.tellg();
  getline(fin, trash);
  //convert to lowercase
  trash = UpperToLower(trash);
  loc = trash.find(elemline);
  if(loc == std::string::npos){
    std::cerr << "CRUNCH ASCII I/O: ELEMENT HEADER NOT IN EXPECTED LOCATION -- "
	      << "remove comments from file (i.e. lines with !)" << std::endl;
    return(1);
  }
  elemloc = fin.tellg();
  //count # of elements
  trash.clear();
  nelem[TET] = 0;
  nelem[PYRAMID] = 0;
  nelem[PRISM] = 0;
  nelem[HEX] = 0;
  numvolelem = 0;
  while(1){
    numvolelem++;
    trash.clear();
    if(fin.peek() >= '0' && fin.peek() <= '9'){
      fin >> elemtype;
      switch (elemtype){
      case CTet:
	nelem[TET]++;
	break;
      case CHex:
	nelem[HEX]++;
	break;
      case CPrism:
	nelem[PRISM]++;
	break;
      case CPyramid:
	nelem[PYRAMID]++;
	break;
      default:
	std::cerr << "ELEM type undefined " << elemtype << std::endl;
      }
    }
    else{
      break;
    }
    getline(fin,trash);
  }
  numvolelem--;
  //stream likely to be in bad state, clear it
  fin.clear();
  lnelem = numvolelem + nelem[TRI] + nelem[QUAD];
    
  std::cout << "CRUNCH ASCII I/O: Number of nodes " << nnode << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Tris " << nelem[TRI] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Quads " << nelem[QUAD] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Tets " << nelem[TET] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Pyramids " << nelem[PYRAMID] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Prisms " << nelem[PRISM] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Hexes " << nelem[HEX] << std::endl;

  //call init routine for mesh
  MemInitMesh();
  std::cout << "CRUNCH ASCII I/O: Memory initialized" << std::endl;

  std::cout << "CRUNCH ASCII I/O: Reading coordinates" << std::endl;

  //read in nodes
  fin.seekg(nodesloc);
  for(i = 0; i < nnode; i++){
    for(j = 0; j < 3; j++){
      fin >> xyz[3*i + j];
    } 
  }
  std::cout << "CRUNCH ASCII I/O: Reading boundary faces" << std::endl;


  //NOTE: CRUNCH nodes are 1 based so they need to be shifted
  //to be zero based... ce la vi
  
  //read in Tris and quads, as well as face tags
  fin.seekg(bfaceloc);
  indx1 = indx2 = 0;
  for(i = 0; i < numbfaces; i++){
    fin >> tag >> npts;
    if(tag < 0){
      std::cout << "WARNING!! negative tags found.. DO NOT USE" << std::endl;
    }
    if(tag == 0){
      tag_shift = 1;
    }
    switch(npts){
      //Triangle
    case 3:
      for(j = 0; j < 3; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[TRI], TRI, 0);
      tempe = new Triangle<Type>;
      tempe->Init(nodes);
      tempe->SetFactag(abs(tag));
      elementList.push_back(tempe);
      indx1++;
      break;
      //Quad
    case 4:
      for(j = 0; j < 4; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[QUAD], QUAD, 0);
      tempe = new Quadrilateral<Type>;
      tempe->Init(nodes);
      tempe->SetFactag(abs(tag));
      elementList.push_back(tempe);
      indx2++;
      break;
      //error
    default:
      std::cerr << "Boundary element type unknown with " << npts << " nodes" << std::endl;
    }
  }
  
  //flip all tags for boundaries to negative
  //and if a zero tag was used shift all tags
  //away -- zero is used for (no tag) condition
  if(tag_shift == 1){
    std::cout << "\n\n";
    std::cout << "CRUNCH ASCII I/O: WARNING!!! ZERO TAGS RESERVED... DO NOT USE" << std::endl;
    std::cout << "CRUNCH ASCII I/O: WARNING!!! ATTEMPTING TO SHIFT TAG IDS BY ONE" << std::endl;
    std::cout << "CRUNCH ASCII I/O: WARNING!!! TAGS WILL BE WRITTEN OUT SHIFTED" << std::endl;
    std::cout << "\n\n";
    for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
      Element<Type>& element = **it;
      Int etype = element.GetType();
      if(etype == TRI || etype == QUAD){
	Int factag = element.GetFactag();
	element.SetFactag(factag + 1);
      }
    }
  }
  std::cout << "CRUNCH ASCII I/O: Reading volume mesh" << std::endl;

  //read in volume elements
  fin.seekg(elemloc);
  indx1 = indx2 = indx3 = indx4 = 0;

  for(i = 0; i < numvolelem; i++){
    fin >> elemtype >> npts;
    switch (elemtype){
    case CTet:
      for(j = 0; j < 4; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[TET], TET, 0);
      tempe = new Tetrahedron<Type>;
      tempe->Init(nodes);
      elementList.push_back(tempe);
      indx1++;
      break;
    case CHex:
      for(j = 0; j < 8; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[HEX], HEX, 0);
      tempe = new Hexahedron<Type>;
      tempe->Init(nodes);
      elementList.push_back(tempe);
      indx2++;
      break;
    case CPyramid:
      for(j = 0; j < 5; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[PYRAMID], PYRAMID, 0);
      tempe = new Pyramid<Type>;
      tempe->Init(nodes);
      elementList.push_back(tempe);
      indx3++;
      break;
    case CPrism:
      for(j = 0; j < 6; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, mnode[PRISM], PRISM, 0);
      tempe = new Prism<Type>;
      tempe->Init(nodes);
      elementList.push_back(tempe);
      indx4++;
      break;
    default:
      std::cerr << "Volume element type unknown -- CRUNCH type " << elemtype << std::endl;
    }
  }

  std::cout << "CRUNCH ASCII I/O: File read successful!!" << std::endl;

  fin.close();
  return(0);
}


template <class Type>
Int Mesh<Type>::WriteXMLVTK_Binary(std::string casename, std::vector<SolutionField<Type>*>& fields)
{
  std::ofstream fout;
  int i, j, e;
  std::string filename = casename + ".vtu";

  //appended data binary blob
  std::vector<unsigned char> appendedData;
  //offsets for each blob to jump into the binary blob
  uint32_t appendedDataAccumulate = 0;
  //data block sizes
  std::vector<uint32_t> appendedDataSizes;
  
  std::cout << "VTK XML BINARY I/O: Writing solution file --> " << filename << std::endl;

  fout.open(filename.c_str());

  //Overall file structure looks like:
  //  <VTKFile type=" UnstructuredGrid" ...> 
  // <UnstructuredGrid> 
  // <Piece NumberOfPoints="#" NumberOfCells="#"> 
  // <PointData Scalars=" Temperature" Vectors=" Velocity"> 
  // <DataArray Name=" Velocity" .../> 
  // <DataArray Name=" Temperature" .../> 
  // <DataArray Name=" Pressure" .../> 
  // </ PointData> 
  // <CellData>...</ CellData> 
  // <Points>...</ Points> 
  //  <Cells>...</ Cells> 
  // </ Piece> 
  // </ UnstructuredGrid> 
  // </ VTKFile> 
  //write the header
  fout << "<?xml version=\"1.0\"?>" << std::endl;
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << nnode << "\" NumberOfCells=\"" << lnelem << "\">" << std::endl;
  fout << "<Points>" << std::endl;
  //write coordinates
  uint32_t datasize = 3*nnode*sizeof(Real);
  //This switch enables/disables base64 inline encoding -- very SLOOOOW
#if 0 
  fout << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << appendedDataAccumulate << "\"/>" << std::endl;
  //resize to totalBytes/sizeof(unsigned char) OR totalBytes/1;
  union doubleUnion_t
  {
    double d;
    char bytes[sizeof(double)];
  };
  doubleUnion_t doubleUnion1;
  appendedDataSizes.push_back(datasize);
  appendedData.resize(appendedDataAccumulate+datasize);
  //iterate of the data one byte at a time
  int byteCounter = 0;
  for(int i = 0; i < 3*nnode; ++i){
    doubleUnion1.d = xyz[i];
    for(int j = 0; j < sizeof(double); ++j){
      appendedData[appendedDataAccumulate + byteCounter + j] = doubleUnion1.bytes[j];
    }
    byteCounter += sizeof(double);
  }
  appendedDataAccumulate += datasize;
#else
  fout << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\" >" << std::endl;
  for(Int i = 0; i < 3*nnode; ++i){
    fout << xyz[i] << " ";
  }
  fout << "</DataArray>" << std::endl;
#endif
  // NOTE: Binary Inline vtu Format
  // The XML structure is like in the ASCII vtu format. But the data is not stored the same
  // way as in the binary vtk format, because the XML language forbids special characters.
  // So the binary data has to be base64-encoded. Additionally, there is a header prepended
  // to the data; it is a 32 bit integer containing the data length (in bytes). This header
  // is encoded separately. So in pseudocode, the data output would look like this
  // (always without spaces or line breaks!):
  //int32 = length( data );
  //output = base64-encode( int32 ) + base64-encode( data );
  // base 64 encode the size of data and then the data itself -
  // the "binary" keyword is really inline base64 which is garbage
  // output = base64-encode( int32 ) + base64-encode( data );
  //fout << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\" \>" << std::endl;
  //std::string datasize64 = base64_encode((const unsigned char*)&datasize, sizeof(uint32_t));
  //std::string data64 = base64_encode((const unsigned char*)xyz, datasize);
  //for(Int i = 0; i < nnode*3; ++i){
  //  fout << datasize64 << data64;
  //}
  //fout << "</DataArray>" << std::endl;

  fout << "</Points>" << std::endl;

  fout << "<Cells>" << std::endl;
 
  // Write out element connectivities
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    for(j = 0; j < nnodes; j++){
      //fout.write((char*)&nodes[j], sizeof(Int));
      fout << nodes[j] << " ";
    }
  }
  fout << "</DataArray>" << std::endl;

  // Write out element offsets into connectivity array
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  Int offset = 0;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    offset += nnodes;
    //fout.write((char*)&offset, sizeof(Int));
    fout << offset << " ";
  }
  fout << "\n</DataArray>" << std::endl;

    // Write out element types
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  Int id;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    e = element.GetType();
    switch (e) {  
    case TRI :
      id = 5;
      //fout.write((char*)&id, sizeof(Int));
      fout << "5" << std::endl;
      break;
    case QUAD :
      id = 9;
      //fout.write((char*)&id, sizeof(Int));
      fout << "9" << std::endl;
      break;
    case TET :
      id = 10;
      //fout.write((char*)&id, sizeof(Int));
      fout << "10" << std::endl;
      break;
    case PYRAMID :
      id = 14;
      //fout.write((char*)&id, sizeof(Int));
      fout << "14" << std::endl;
      break;
    case PRISM :
      id = 13;
      //fout.write((char*)&id, sizeof(Int));
      fout << "13" << std::endl;
      break;
    case HEX :
      id = 12;
      //fout.write((char*)&id, sizeof(Int));
      fout << "12" << std::endl;
      break;
    default :
      std::cerr << "Type not defined WriteXMLVTK_BINARY() -- type " << e << std::endl;
    }
  }
  fout << "\n</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  
  
  fout << "<CellData>" << std::endl;
  fout << "<DataArray type=\"Float64\" Name=\"BoundaryId\" format=\"ascii\">" << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    e = element.GetFactag();
    //fout.write((char*)&rint, sizeof(Int));
    fout << e << " ";
  }
  fout << "</DataArray>" << std::endl;
  fout << "</CellData>" << std::endl;
  if(fields.size() != 0){
    fout << "<PointData>" << std::endl;
    for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){ 
      SolutionField<Type>& field = **it;
      Real* q = field.GetData(FIELDS::STATE_NONE);
      Int neqn = field.GetNdof();
      for(j = 0; j < neqn; j++){
	if(field.DofIsVector(j)){
	  fout << "<DataArray type=\"Float64\" Name=\"" << field.GetDofName(j) << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	  for(i = 0; i < nnode; i++){
	    //we check for nan's here for a few reasons... mostly because paraview chokes if we dump them out but with
	    //a really oddball array not long enough error which is not useful. So we do something clever and dump a huge actual number
	    //instead b/c at least that allows diagnostics
	    if(isnan(real(q[i*neqn + j + 0])) || isnan(real(q[i*neqn + j +1])) || isnan(real(q[i*neqn + j + 2]))){
	      fout << std::numeric_limits<double>::max() << " " << std::numeric_limits<double>::max() << " " <<std::numeric_limits<double>::max() << std::endl;
	    }
	    else{
	      //fout.write((char*)&q[i*neqn + j], sizeof(Real)*3);
	      fout << q[i*neqn + j + 0] << " " << q[i*neqn + j + 1] << " " << q[i*neqn + j + 2] << " " ;
	    }
	  }
	  j+=2;
	  fout << "</DataArray>";
	}
	else{
	  fout << "<DataArray type=\"Float64\" Name=\"" << field.GetDofName(j) << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for(i = 0; i < nnode; i++){
	    //fout.write((char*)&q[i*neqn + j], sizeof(Real));
	    fout << q[i*neqn + j] << " ";
	  }
	  fout << "</DataArray>";
	}
      }
    }
    fout << "</PointData>";
  }

  //write file closing
  fout << "</Piece>" << std::endl;
  fout << "</UnstructuredGrid>" << std::endl;


#if 0
  //write the appended data block
  //The format is
  // _NNNN<data>NNNN<data>NNNN<data>
  // ^         ^         ^
  // 1         2         3
  //where each "NNNN" is an unsigned 32-bit integer, and <data> consists of
  //a number of bytes equal to the preceding NNNN value.  The corresponding
  //DataArray elements must have format="appended" and offset attributes
  //equal to the following:
  //1.) offset="0"
  //2.) offset="(4+NNNN1)"
  //3.) offset="(4+NNNN1+4+NNNN2)"
  fout << "<AppendedData encoding=\"raw\">" << std::endl;
  fout << "_";
  fout.write((char*)&appendedDataSizes[0], sizeof(uint32_t));
  fout.write((const char*)appendedData.data(), appendedData.size());
  fout << "</AppendedData>" << std::endl;
#endif
  fout << "</VTKFile>" << std::endl;
  
  fout.close();

  std::cout << "VTK XML BINARY I/O: File write successful!!" << std::endl;

  return (0);

}

template <class Type>
Int Mesh<Type>::WriteVTK_Binary(std::string casename, std::vector<SolutionField<Type>*>& fields)
{
  std::ofstream fout;
  int i, j, e;
  std::string filename = casename + ".vtk";

  std::cout << "VTK BINARY I/O: Writing solution file --> " << filename << std::endl;

  fout.open(filename.c_str());

  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << "Solver solution data" << std::endl;
  fout << "BINARY" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  fout << "POINTS " << nnode << " double" << std::endl;
  
  //write coordinates
  ReverseEndianArrayReal(xyz, nnode*3);
  fout.write((char*)xyz, sizeof(Real)*nnode*3);
  ReverseEndianArrayReal(xyz, nnode*3);
  fout << "\n";
  //fout << xyz[i*3 + 0] << " " << xyz[i*3 + 1] << " " << xyz[i*3 + 2] << std::endl;
  
  j = 0;
  for(e = TRI; e <= HEX; e++){
    j += nelem[e]*mnode[e];
  }
  j += lnelem;
  fout << "CELLS " << lnelem << " " << j << std::endl; 
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    Int rint = ReverseEndianInt(nnodes);
    fout.write((char*)&rint, sizeof(Int));
    for(j = 0; j < nnodes; j++){
      rint = ReverseEndianInt(nodes[j]);
      fout.write((char*)&rint, sizeof(Int));
    }
  }
  fout << "\n";

  fout << "CELL_TYPES " << lnelem << std::endl;
  Int id;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    e = element.GetType();
    switch (e) {  
    case TRI :
      id = 5;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "5" << std::endl;
      break;
    case QUAD :
      id = 9;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "9" << std::endl;
      break;
    case TET :
      id = 10;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "10" << std::endl;
      break;
    case PYRAMID :
      id = 14;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "14" << std::endl;
      break;
    case PRISM :
      id = 13;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "13" << std::endl;
      break;
    case HEX :
      id = 12;
      id = ReverseEndianInt(id);
      fout.write((char*)&id, sizeof(Int));
      //fout << "12" << std::endl;
      break;
    default :
      std::cerr << "Type not defined WriteVTK_BINARY() -- type " << e << std::endl;
    }
  }
  
  fout << "\n";
  
  
  fout << "CELL_DATA " << lnelem << std::endl;
  fout << "SCALARS BoundaryIds int " << 1 << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for(typename std::vector<Element<Type>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Type>& element = **it;
    e = element.GetFactag();
    Int rint = ReverseEndianInt(e);
    fout.write((char*)&rint, sizeof(Int));
  }
  fout << "\n";
  if(fields.size() != 0){
    fout << "POINT_DATA " << nnode << std::endl;
  }
  for(typename std::vector<SolutionField<Type>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Type>& field = **it;
    Real* q = field.GetData(FIELDS::STATE_NONE);
    Int neqn = field.GetNdof();
    for(j = 0; j < neqn; j++){
      if(field.DofIsVector(j)){
	fout << "VECTORS " << field.GetDofName(j) << " double " << std::endl;
	for(i = 0; i < nnode; i++){
	  ReverseEndianArrayReal(&q[i*neqn + j], 3);
	  fout.write((char*)&q[i*neqn + j], sizeof(Real)*3);
	  ReverseEndianArrayReal(&q[i*neqn + j], 3);
	}
	j+=2;
	fout << "\n";
      }
      else{
	fout << "SCALARS " << field.GetDofName(j) << " double " << 1 << std::endl;
	fout << "LOOKUP_TABLE default" << std::endl;
	for(i = 0; i < nnode; i++){
	  Real temp = ReverseEndianReal(q[i*neqn + j]);
	  fout.write((char*)&temp, sizeof(Real));
	}
	fout << "\n";
      }
    }
  }



  fout.close();

  std::cout << "VTK BINARY I/O: File write successful!!" << std::endl;


  return (0);
}

