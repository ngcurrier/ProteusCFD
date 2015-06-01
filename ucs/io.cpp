#include "io.h"
#include "etypes.h"
#include "endian_util.h"
#include "strings_util.h"
#include "tinyxml.h"
#include "parallel.h"
#include "h5layer.h"
#include "mesh.h"
#include "dataInfo.h"
#include "solutionField.h"

void TranslateWinding(Int* nodes, Int translation[6][8], Int num_nodes, Int etype, Int to_other_format){
  Int tempnode[num_nodes];
  Int k;

  memcpy(tempnode,nodes,num_nodes*sizeof(Int));
  
  if(to_other_format){
    for (k = 0; k < num_nodes; k++){
      nodes[k] = tempnode[translation[etype][k]];
    }
  }
  else{
    for (k = 0; k < num_nodes; k++){
      // translation[k] gives our local node number
      // k = other format local node number
      nodes[translation[etype][k]] = tempnode[k];
    }
  }
}

Int ReadUGRID_Ascii(Mesh<Real> &m, std::string filename)
{
  std::ifstream fin;
  Int i, j;
  Bool tag_shift = 0;
  Int nodes[8];
  Element<Real>* tempe = NULL;

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

  fin >> m.nnode >> m.nelem[TRI] >> m.nelem[QUAD] >> m.nelem[TET];
  fin >> m.nelem[PYRAMID] >> m.nelem[PRISM] >> m.nelem[HEX];

  std::cout << "UGRID ASCII I/O: Number of nodes " << m.nnode << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Tris " << m.nelem[TRI] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Quads " << m.nelem[QUAD] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Tets " << m.nelem[TET] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Pyramids " << m.nelem[PYRAMID] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Prisms " << m.nelem[PRISM] << std::endl;
  std::cout << "UGRID ASCII I/O: Number of Hexes " << m.nelem[HEX] << std::endl;

  m.lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    m.lnelem += m.nelem[e];
  }

  //call init routine for mesh
  m.MemInitMesh();
  std::cout << "UGRID ASCII I/O: Memory initialized" << std::endl;

  std::cout << "UGRID ASCII I/O: Reading coordinates" << std::endl;

  //read in nodes
  for(i = 0; i < m.nnode; i++){
    for(j = 0; j < 3; j++){
      fin >> m.xyz[3*i + j];
    } 
  }
  std::cout << "UGRID ASCII I/O: Reading boundary faces" << std::endl;


  //NOTE: UGRID nodes are 1 based so they need to be shifted
  //to be zero based... ce la vi
  
  //read in Tris
  for(i = 0; i < m.nelem[TRI]; i++){
    for(j = 0; j < 3; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[TRI], TRI, 0);
    tempe = new Triangle<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }
  //read in Quads
  //UGRID winds these backwards
  for(i = 0; i < m.nelem[QUAD]; i++){
    for(j = 0; j < 4; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[QUAD], QUAD, 0);
    tempe = new Quadrilateral<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }
  
  //read in boundary tags for tris
  m.nfactags = 0;
  Int ielem = 0;
  for(i = 0; i < m.nelem[TRI]; i++){
    Int temp;
    fin >> temp;
    if(temp < 0){
      std::cout << "WARNING!! negative tags found.. DO NOT USE" << std::endl;
    }
    if(temp == 0){
      tag_shift = 1;
    }
    m.elementList[ielem]->SetFactag(abs(temp));
    m.nfactags = MAX(m.nfactags, m.elementList[ielem]->GetFactag());
    ielem++;
  }
  //read in boundary tags for quads
  for(i = 0; i < m.nelem[QUAD]; i++){
    Int temp;
    fin >> temp;
    if(temp < 0){
      std::cout << "WARNING!! negative tags found.. DO NOT USE" << std::endl;
    }
    if(temp == 0){
      tag_shift = 1;
    }
    m.elementList[ielem]->SetFactag(abs(temp));
    m.nfactags = MAX(m.nfactags, m.elementList[ielem]->GetFactag());
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
    for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
      Element<Real>& element = **it;
      Int etype = element.GetType();
      if(etype == TRI || etype == QUAD){
	Int factag = element.GetFactag();
	element.SetFactag(factag + 1);
      }
    }
  }
  std::cout << "UGRID ASCII I/O: Reading volume mesh" << std::endl;

  //read in Tets
  for(i = 0; i < m.nelem[TET]; i++){
    for(j = 0; j < 4; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[TET], TET, 0);
    tempe = new Tetrahedron<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }
  //read in Pyramids
  //UGRID winds these weird... we don't use their winding
  for(i = 0; i < m.nelem[PYRAMID]; i++){
    for(j = 0; j < 5; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[PYRAMID], PYRAMID, 0);
    tempe = new Pyramid<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }
  //read in Prisms
  for(i = 0; i < m.nelem[PRISM]; i++){
    for(j = 0; j < 6; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[PRISM], PRISM, 0);
    tempe = new Prism<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }
  //read in Hexes
  for(i = 0; i < m.nelem[HEX]; i++){
    for(j = 0; j < 8; j++){
      fin >> nodes[j];
      nodes[j]--;
    }
    TranslateWinding(nodes, translation, m.mnode[HEX], HEX, 0);
    tempe = new Hexahedron<Real>;
    tempe->Init(nodes);
    m.elementList.push_back(tempe);
  }

  std::cout << "UGRID ASCII I/O: File read successful!!" << std::endl;

  fin.close();
  return(0);
}

Int ReadCRUNCH_Ascii(Mesh<Real> &m, std::string filename)
{
  Element<Real>* tempe = NULL;
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
  fin >> m.nfactags;
  getline(fin, trash);
  //skip boundary table lines
  for(i = 0; i < m.nfactags; i++){
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
  fin >> m.nnode;
  getline(fin, trash);
  nodesloc = fin.tellg();
  //burn through node lines... we'll read them later
  for(i = 0; i < m.nnode; i++){
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
  m.nelem[TRI] = 0;
  m.nelem[QUAD] = 0;
  for(i = 0; i < numbfaces; i++){
    fin >> temp >> npts;
    getline(fin, trash);
    if(npts == 3){
      m.nelem[TRI]++;
    }
    else if(npts == 4){
      m.nelem[QUAD]++;
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
  m.nelem[TET] = 0;
  m.nelem[PYRAMID] = 0;
  m.nelem[PRISM] = 0;
  m.nelem[HEX] = 0;
  numvolelem = 0;
  while(1){
    numvolelem++;
    trash.clear();
    if(fin.peek() >= '0' && fin.peek() <= '9'){
      fin >> elemtype;
      switch (elemtype){
      case CTet:
	m.nelem[TET]++;
	break;
      case CHex:
	m.nelem[HEX]++;
	break;
      case CPrism:
	m.nelem[PRISM]++;
	break;
      case CPyramid:
	m.nelem[PYRAMID]++;
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
  m.lnelem = numvolelem + m.nelem[TRI] + m.nelem[QUAD];
    
  std::cout << "CRUNCH ASCII I/O: Number of nodes " << m.nnode << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Tris " << m.nelem[TRI] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Quads " << m.nelem[QUAD] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Tets " << m.nelem[TET] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Pyramids " << m.nelem[PYRAMID] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Prisms " << m.nelem[PRISM] << std::endl;
  std::cout << "CRUNCH ASCII I/O: Number of Hexes " << m.nelem[HEX] << std::endl;

  //call init routine for mesh
  m.MemInitMesh();
  std::cout << "CRUNCH ASCII I/O: Memory initialized" << std::endl;

  std::cout << "CRUNCH ASCII I/O: Reading coordinates" << std::endl;

  //read in nodes
  fin.seekg(nodesloc);
  for(i = 0; i < m.nnode; i++){
    for(j = 0; j < 3; j++){
      fin >> m.xyz[3*i + j];
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
      TranslateWinding(nodes, translation, m.mnode[TRI], TRI, 0);
      tempe = new Triangle<Real>;
      tempe->Init(nodes);
      tempe->SetFactag(abs(tag));
      m.elementList.push_back(tempe);
      indx1++;
      break;
      //Quad
    case 4:
      for(j = 0; j < 4; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, m.mnode[QUAD], QUAD, 0);
      tempe = new Quadrilateral<Real>;
      tempe->Init(nodes);
      tempe->SetFactag(abs(tag));
      m.elementList.push_back(tempe);
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
    for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
      Element<Real>& element = **it;
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
      TranslateWinding(nodes, translation, m.mnode[TET], TET, 0);
      tempe = new Tetrahedron<Real>;
      tempe->Init(nodes);
      m.elementList.push_back(tempe);
      indx1++;
      break;
    case CHex:
      for(j = 0; j < 8; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, m.mnode[HEX], HEX, 0);
      tempe = new Hexahedron<Real>;
      tempe->Init(nodes);
      m.elementList.push_back(tempe);
      indx2++;
      break;
    case CPyramid:
      for(j = 0; j < 5; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, m.mnode[PYRAMID], PYRAMID, 0);
      tempe = new Pyramid<Real>;
      tempe->Init(nodes);
      m.elementList.push_back(tempe);
      indx3++;
      break;
    case CPrism:
      for(j = 0; j < 6; j++){
	fin >> nodes[j];
	nodes[j]--;
      }
      TranslateWinding(nodes, translation, m.mnode[PRISM], PRISM, 0);
      tempe = new Prism<Real>;
      tempe->Init(nodes);
      m.elementList.push_back(tempe);
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

Int ReadSU2_Ascii(Mesh<Real>& m, std::string filename)
{
  Element<Real>* tempe = NULL;
  std::ifstream fin;
  std::string line;
  Int ndim = 0;   //number of dimensions
  Int nelem = 0;  //number of elements
  Int nnodes = 0; //number of points
  Int nbound = 0; //number of boundaries
    
  Int elemCounter = 0; //counts number of elements read
  Int nodeCount = 0; //counts number of nodes read
  Int boundElemCounter = 0; // counts boundary elements read
  Int boundCounter = 0; // counts boundary sections read
  Int elemType;
  Int nbelemTmp = 0; // temporary counter for each boundary section

  //initialize element type counters
  m.nelem[TRI] = m.nelem[QUAD] = m.nelem[TET] = m.nelem[PYRAMID] = m.nelem[PRISM] = m.nelem[HEX] = 0;

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
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NDIME=\" marker where expected";
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
	  ss >> nelem;
	  std::cout << "SU2 ASCII I/O: Reading " << nelem << " elements" << std::endl;
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
	  tempe = new Quadrilateral<Real>;
	  tempe->Init(nodes);
	  m.elementList.push_back(tempe);
	  m.nelem[TET]++;
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
	  tempe = new Hexahedron<Real>;
	  tempe->Init(nodes);
	  m.elementList.push_back(tempe);
	  m.nelem[HEX]++;
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
	  tempe = new Prism<Real>;
	  tempe->Init(nodes);
	  m.elementList.push_back(tempe);
	  m.nelem[PRISM]++;
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
	  tempe = new Pyramid<Real>;
	  tempe->Init(nodes);
	  m.elementList.push_back(tempe);
	  m.nelem[PYRAMID]++;
	}
	else{
	  std::cerr << "SU2 ASCII I/O: Cannot find element of type " << elemType << std::endl;
	  return(-1);
	}
	elemCounter++;
	if(elemCounter >= nelem){
	  state = stateReadNpoints;
	}
	break;
      case stateReadNpoints:
	loc = line.find(npoinKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NPOIN=\" marker where expected";
	  return(1);
	}
	else{
	  ss.str(line.substr(loc+npoinKey.size()));
	  ss >> nnodes;
	  std::cout << "SU2 ASCII I/O: Reading " << nnodes << " nodes" << std::endl;
	  state = stateReadPoints;
	  m.nnode = nnodes;
	  m.MemInitMesh();
	  std::cout << "SU2 ASCII I/O: Memory initialized" << std::endl;
	}
	break;
      case stateReadPoints:
	ss.str(line);
	{
	  Real xyz[3];
	  ss >> xyz[0];
	  ss >> xyz[1];
	  ss >> xyz[2];
	  Int id;
	  ss >> id;
	  m.xyz[nodeCount*3 + 0] = xyz[0];
	  m.xyz[nodeCount*3 + 1] = xyz[1];
	  m.xyz[nodeCount*3 + 2] = xyz[2];
	  nodeCount++;
	}
	if(nodeCount >= nnodes){
	  state = stateReadNmark;
	}
	break;
      case stateReadNmark:
	loc = line.find(nmarkKey);
	if(loc == std::string::npos){
	  std::cerr << "SU2 ASCII I/O: File read did not find \"NMARK=\" marker where expected";
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
	  std::cerr << "SU2 ASCII I/O: File read did not find \"MARKER_TAG=\" where expected";
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
	  std::cerr << "SU2 ASCII I/O: File read did not find \"MARKER_ELEMS=\" where expected";
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
	  std::cerr << "SU2 ASCII I/O: Only 3d elements expected, found line";
	  return (-1);
	}
	else if(elemType == SU2_TRI){
	  Int nodes[3];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  Int id;
	  ss >> id;
	  tempe = new Triangle<Real>;
	  tempe->Init(nodes);
	  tempe->SetFactag(boundCounter+1);
	  m.elementList.push_back(tempe);
	  m.nelem[TRI]++;
	}
	else if(elemType == SU2_QUAD){
	  Int nodes[4];
	  ss >> nodes[0];
	  ss >> nodes[1];
	  ss >> nodes[2];
	  ss >> nodes[3];
	  Int id;
	  ss >> id;
	  tempe = new Quadrilateral<Real>;
	  tempe->Init(nodes);
	  tempe->SetFactag(boundCounter+1);
	  m.elementList.push_back(tempe);
	  m.nelem[QUAD]++;
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

  m.nfactags = nbound;
  m.lnelem = 0;
  for(Int e = TRI; e <= HEX; e++){
    m.lnelem += m.nelem[e];
  }
  std::cout << "SU2 ASCII I/O: Number of nodes " << m.nnode << std::endl;
  std::cout << "SU2 ASCII I/O: Number of elements " << m.lnelem << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Tris " << m.nelem[TRI] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Quads " << m.nelem[QUAD] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Tets " << m.nelem[TET] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Pyramids " << m.nelem[PYRAMID] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Prisms " << m.nelem[PRISM] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of Hexes " << m.nelem[HEX] << std::endl;
  std::cout << "SU2 ASCII I/O: Number of tagged boundaries " << m.nfactags << std::endl;

  fin.close();
  return(0);
}

Int WriteCRUNCH_Ascii(Mesh<Real> &m, std::string casename)
{
  Int i,j;
  Int etype;    
  Int tempnodes[8];

  std::ofstream fout;
  std::stringstream ss (std::stringstream::in | std::stringstream::out);

  std::string filename = casename + ".CRUNCH";
  
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
  fout << "Boundary Table " << m.nfactags << std::endl;
  for(i = 0; i < m.nfactags; i++){
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
  
  fout << "NODES " << m.nnode << std::endl;

  //write nodes
  fout.setf(std::ios::scientific);
  fout.precision(15);
  for(i = 0; i < m.nnode; i++){
    fout << m.xyz[i*3 + 0] << " " << m.xyz[i*3 + 1] << " " << m.xyz[i*3 + 2] << std::endl;
  }
  
  //NOTE: CRUNCH nodes are 1 based so they need to be shifted
  //from zero base used internally

  //write boundary faces
  fout << "BOUNDARY FACES " << m.nelem[TRI]+m.nelem[QUAD] << std::endl;
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
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
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
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

Int WriteVTK_Ascii(Mesh<Real> &m, std::string casename, std::vector<SolutionField<Real>*>& fields)
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
  fout << "POINTS " << m.nnode << " double" << std::endl;
  for(Int i = 0; i < m.nnode; i++){
    fout << m.xyz[i*3 + 0] << " " << m.xyz[i*3 + 1] << " " << m.xyz[i*3 + 2] << std::endl;
  }
  j = 0;
  for(Int e = TRI; e <= HEX; e++){
    j += m.nelem[e]*m.mnode[e];
  }
  j += m.lnelem;
  fout << "CELLS " << m.lnelem << " " << j << std::endl; 
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
    Int* nodes = NULL;
    Int nnodes = element.GetNodes(&nodes);
    fout << nnodes << " ";
    for(Int i = 0; i < nnodes; i++){
      fout << nodes[i] << " ";
    }
    fout << "\n";
  }

  fout << "CELL_TYPES " << m.lnelem << std::endl;
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
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
  
  fout << "CELL_DATA " << m.lnelem << std::endl;
  fout << "SCALARS BoundaryIds int " << 1 << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
    Int e = element.GetFactag();
    fout << e << "\n" ;
  }
  if(fields.size() != 0){
    fout << "POINT_DATA " << m.nnode << std::endl;
  }
  for(std::vector<SolutionField<Real>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Real>& field = **it;
    Int neqn = field.GetNdof();
    Real* q = field.GetData(FIELDS::STATE_NONE);
    for(j = 0; j < neqn; j++){
      if(field.DofIsVector(j)){
	fout << "VECTORS " << field.GetDofName(j) << " double " << std::endl;
	for(Int i = 0; i < m.nnode; i++){
	  fout << q[i*neqn + j] << " ";
	  fout << q[i*neqn + j+1] << " ";
	  fout << q[i*neqn + j+2] << std::endl;;
	}
	j+=2;
      }
      else{
	fout << "SCALARS " << field.GetDofName(j) << " double " << 1 << std::endl;
	fout << "LOOKUP_TABLE default" << std::endl;
	for(Int i = 0; i < m.nnode; i++){
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

Int WriteVTK_Binary(Mesh<Real> &m, std::string casename, std::vector<SolutionField<Real>*>& fields)
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
  fout << "POINTS " << m.nnode << " double" << std::endl;
  
  //write coordinates
  ReverseEndianArrayReal(m.xyz, m.nnode*3);
  fout.write((char*)m.xyz, sizeof(Real)*m.nnode*3);
  ReverseEndianArrayReal(m.xyz, m.nnode*3);
  fout << "\n";
  //fout << m.xyz[i*3 + 0] << " " << m.xyz[i*3 + 1] << " " << m.xyz[i*3 + 2] << std::endl;
  
  j = 0;
  for(e = TRI; e <= HEX; e++){
    j += m.nelem[e]*m.mnode[e];
  }
  j += m.lnelem;
  fout << "CELLS " << m.lnelem << " " << j << std::endl; 
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
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

  fout << "CELL_TYPES " << m.lnelem << std::endl;
  Int id;
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
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
  
  
  fout << "CELL_DATA " << m.lnelem << std::endl;
  fout << "SCALARS BoundaryIds int " << 1 << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for(std::vector<Element<Real>*>::iterator it = m.elementList.begin(); it != m.elementList.end(); ++it){
    Element<Real>& element = **it;
    e = element.GetFactag();
    Int rint = ReverseEndianInt(e);
    fout.write((char*)&rint, sizeof(Int));
  }
  fout << "\n";
  if(fields.size() != 0){
    fout << "POINT_DATA " << m.nnode << std::endl;
  }
  for(std::vector<SolutionField<Real>*>::iterator it = fields.begin(); it != fields.end(); ++it){
    SolutionField<Real>& field = **it;
    Real* q = field.GetData(FIELDS::STATE_NONE);
    Int neqn = field.GetNdof();
    for(j = 0; j < neqn; j++){
      if(field.DofIsVector(j)){
	fout << "VECTORS " << field.GetDofName(j) << " double " << std::endl;
	for(i = 0; i < m.nnode; i++){
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
	for(i = 0; i < m.nnode; i++){
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

Int WriteGridXDMF(Mesh<Real> &m, PObj<Real> &p, std::string filebase, std::string meshname)
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
  
  nnode[0] = m.nnode+m.gnode;

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
  HDF_WriteArray(fileId, directory, dataname, m.xyz, (m.nnode+m.gnode)*3);
						     
  HDF_CloseFile(fileId);

  delete [] elements;

  return err;
}
