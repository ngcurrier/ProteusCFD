#include "mesh.h"
#include "general.h"
#include "h5layer.h"
#include "mem_util.h"
#include "macros.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <metis.h>
#include <vector>

using namespace std;

void METISpartition(Mesh<Real>* m, Int* part, Int np, Int* nodesPart);
void Coordpartition(Mesh<Real>* m, Int* part, Int np, Int* nodesPart);

int main(int argc, char* argv[]){

  HDF_TurnOffErrorHandling();

  //mesh has not been reordered at this point
  Bool reordered = false;

  Mesh<Real> m;
  string filename;
  string fileextension;
  string casename;
  size_t pos;
  Int np = 0;
  
  //use these as dynamic holders, fast and easy
  vector<Int> vint;
  vector<Real> vreal;

  Int i, j, k;
  Int ierr = 0;

  if(argc != 3){
    cerr << "Invalid arguments" << endl;
    cerr << "<Mesh file> <# of partitions>" << endl;
    return(3);
  }

  np = atoi(argv[2]);
  if(np <= 0){
    cerr << "Number of partitions must be > 0" << endl;
    return(4);
  }

  //input file name parsing -- reader selection
  filename = argv[1];
  pos = filename.rfind(".");
  casename = filename.substr(0,pos);
  fileextension = filename.substr(pos+1, filename.size() - (casename.size() ));
  if(fileextension == "crunch"){
    ierr = m.ReadCRUNCH_Ascii(filename);
  }
  else if(fileextension == "ugrid"){
    ierr = m.ReadUGRID_Ascii(filename);
  }
  else if(fileextension == "su2"){
    ierr = m.ReadSU2_Ascii(filename);
  }
  else if(fileextension == "msh"){
    ierr = m.ReadGMSH_Ascii(filename);
  }
  else{
    cerr << "File extension " << fileextension << " not recognized" << endl;
    return(0);
  }

  if(ierr){
    cerr << "Grid read failure: could not open -- " << filename << endl;
    return (1);
  }

  //grab node coords
  double const * xyz = m.GetNodeCoords();
  
  //build utility maps for METIS
  m.BuildMapsDecomp();

  Int nnode = m.GetNumNodes();
  Int lnelem = m.GetNumLocalElem();
  idx_t* part = new idx_t[nnode];
  Int localNodesCount[np];
  Int localNodesOffset[np];
  Int ghostNodesCount[np];
  Int localElemsCount[np];
  Int splitElemsCount[np];

  for(i = 0; i < np; i++){
    localNodesCount[i] = 0;
    ghostNodesCount[i] = 0;
    localElemsCount[i] = 0;
    splitElemsCount[i] = 0;
  }

  //get partition vector
  //localNodesCount returns correct from this routine
  METISpartition(&m, part, np, localNodesCount);
  //Coordpartition(&m, part, np, localNodesCount);

  //set up local node offsets
  localNodesOffset[0] = 0;
  for(i = 1; i < np; i++){
    localNodesOffset[i] = localNodesOffset[i-1] + localNodesCount[i-1];
  }

  //allocate and load up local nodes global Id
  Int* localNodes = new Int[nnode];
  for(i = 0; i < np; i++){
    Int count = 0;
    for(j = 0; j < nnode; j++){
      if(part[j] == i){
	localNodes[localNodesOffset[i] + count] = j;
	count++;
      }
    }
  }

  //assume the likely number of ghost nodes drops with number of domains
  Int goffset = nnode/np + 5;     
  Int* ghostNodes = new Int[np * goffset];
  Int* ghostDoms = new Int[goffset];
  Int leoffset[np]; 
  Int** localElems = new Int*[np];
  for(i = 0; i < np; i++){
    //add an extra 20% as a buffer for unbalanced split
    leoffset[i] = lnelem/np + lnelem/5;
    localElems[i] = new Int[leoffset[i]];
  }
  //same thing for split elements but with a little wiggle room
  Int seoffset = lnelem/np + 35;
  Int* splitElems = new Int[np * seoffset];

  vint.reserve(nnode/np);
  vreal.reserve(nnode/np*3);
	       
  Bool* elemChecked = new Bool[lnelem];
  for(i = 0; i < lnelem; i++){
    elemChecked[i] = false;
  } 
 
  //count elements that are totally within a single domain
  for(UInt gelem = 0; gelem < m.elementList.size(); gelem++){
    Element<Real>& element = *m.elementList[gelem];
    Int* nodes;
    Int nnodes = element.GetNodes(&nodes);
    Int node = nodes[0];
    Int dom = part[node];
    for(i = 0; i < nnodes; i++){
      if(part[nodes[i]] != dom){
	node = -1;
	break;
      }
    }
    if(node != -1){
      elemChecked[gelem] = true;
      localElems[dom][localElemsCount[dom]] = gelem;
      localElemsCount[dom]++;
      if(localElemsCount[dom] >= leoffset[dom]){
	MemResize(&localElems[dom], leoffset[dom], leoffset[dom]+leoffset[dom]/5);
	leoffset[dom] += leoffset[dom]/5;
      }
    }
  }
  
  //count elements which are unaccounted for and add them to the split 
  //elements list of appropriate domains
  for(i = 0; i < lnelem; i++){
    //if elem was not local to a single domain
    if(elemChecked[i] == false){
      Int inDomain[np];
      for(j = 0; j < np; j++){
	inDomain[j] = 0;
      }
      Element<Real>& element = *m.elementList[i];
      Int* nodes;
      Int nnodes = element.GetNodes(&nodes);
      for(j = 0; j < nnodes; j++){
	Int node = nodes[j];
	Int dom = part[node];
	inDomain[dom] = 1;
      }
      for(j = 0; j < np; j++){
	if(inDomain[j] == 1){
	  splitElems[j*seoffset + splitElemsCount[j]] = i;
	  splitElemsCount[j]++;
	  if(splitElemsCount[j] >= seoffset){
	    cerr << "WARNING: not enough memory for split elem list on domain " << j << endl;
	    cerr << "BUMP seoffset in code -- known shortfall!!!" << endl;
	    cerr << "seoffset value is " << seoffset << endl;
	    return (1);
	  }
	}
      }
    }
  }

  delete [] elemChecked;
	
  //used to make sure duplicates are not put in the ghost nodes list
  //this would wreak all sorts of havoc
  Bool* nodeInList = new Bool[nnode];

  //get ghost nodes list for split elements
  for(i = 0; i < np; i++){
    for(j = 0; j < nnode; j++){
      nodeInList[j] = 0;
    }
    for(j = 0; j < splitElemsCount[i]; j++){
      Int gelem = splitElems[i*seoffset + j];
      Element<Real>& element = *m.elementList[gelem];
      Int* nodes;
      Int nnodes = element.GetNodes(&nodes);
      for(k = 0; k < nnodes; k++){
	Int node = nodes[k];
	if(part[node] != i && nodeInList[node] == 0){
	  nodeInList[node] = 1;
	  ghostNodes[i*goffset + ghostNodesCount[i]] = node;
	  ghostNodesCount[i]++;
	  if(ghostNodesCount[i] >= goffset){
	    cerr << "WARNING: not enough memory for ghost nodes list on domain " << i << endl;
	    cerr << "BUMP goffset in code -- known shortfall!!!" << endl;
	    return (1);
	  }
	}
      }
    }
  }
  delete [] nodeInList;

  //order ghost nodes list so that they are grouped by process
  //i.e. proc 1 ghostnodes, proc2 ghostnodes, etc.
  //makes receiving comms less memory intensive during runtime
  //eat the cost here instead
  Int* tempNodes;
  Int maxGhostNodes = 0;
  for(i = 0; i < np; i++){
    if(ghostNodesCount[i] > maxGhostNodes){
      maxGhostNodes = ghostNodesCount[i];
    }
  }
  tempNodes = new Int[maxGhostNodes];
  for(i = 0; i < np; i++){
    Int count = 0;
    for(j = 0; j < np; j++){
      for(k = 0; k < ghostNodesCount[i]; k++){
	Int node = ghostNodes[i*goffset + k];
	if(part[node] == j){
	  tempNodes[count] = node;
	  count++;
	}
      }
    }
    memcpy(&ghostNodes[i*goffset + 0], tempNodes, sizeof(Int)*count);
    if(count != ghostNodesCount[i]){
      cout << "WARNING -- inconsistent ghost nodes count" << endl;
    }
  }
  delete [] tempNodes;


  /////////////////////////////////////////////////////////////////////////////////
  //
  //      START ---- WRITE DATA FOR SPLIT ELEMENTS/NODES/ETC. TO FILE
  //
  ////////////////////////////////////////////////////////////////////////////////

  ostringstream temposs;
  //open appropriate files
  std::string directoryBase = "";
  hid_t h5out[np];
  for(i = 0; i < np; i++){
    temposs.str("");
    temposs << i;
    std::string h5filename = casename + "." + temposs.str() + ".h5";
    std::cout << "Removing old partition file " + h5filename  + "\n";
    std::string rmline = "rm " + h5filename;
    system(rmline.c_str());
    h5out[i] = HDF_OpenFile(h5filename, 1);
    if(h5out[i] < 0){
      std::cerr << "DECOMP I/O: Could not open file -- " << h5filename << std::endl;
    }
  }
  HDF_TurnOffErrorHandling();
  for(i = 0; i < np; i++){
    HDF_WriteScalar(h5out[i], "/", "Number Of Processors", &np);
    Int reorder = reordered;
    HDF_WriteScalar(h5out[i], "/", "Reordered", &reorder);
    Int tnnode = nnode;
    HDF_WriteScalar(h5out[i], "/", "Global Number Of Nodes", &tnnode);
    Int tlnelem = lnelem;
    HDF_WriteScalar(h5out[i], "/", "Global Number Of Elements", &tlnelem);
    HDF_WriteScalar(h5out[i], "/", "Number Of Factags", &m.nfactags);
    Int scaled = 0;
    HDF_WriteScalar(h5out[i], "/", "Rescaled", &scaled);
  }
  ////////////////////////////////////////////////////////
  //
  //  MESH TOTALS WRITING
  //
  ////////////////////////////////////////////////////////
  
  for(i = 0; i < np; i++){
    directoryBase = "/Mesh/";

    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Local Nodes", &localNodesCount[i]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Ghost Nodes", &ghostNodesCount[i]);

    //write number of each element type including split elements
    Int type_count[6];
    for(j = 0; j < 6; j++){
      type_count[j] = 0;
    }
    for(j = 0; j < localElemsCount[i]; j++){
      Int gelem = localElems[i][j];
      Element<Real>& element = *m.elementList[gelem];
      Int etype = element.GetType();
      type_count[etype]++;
    }
    for(j = 0; j < splitElemsCount[i]; j++){
      Int gelem = splitElems[i*seoffset + j];
      Element<Real>& element = *m.elementList[gelem];
      Int etype = element.GetType();
      type_count[etype]++;
    }
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Triangles", &type_count[TRI]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Quadrilaterals", &type_count[QUAD]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Tetrahedron", &type_count[TET]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Pyramids", &type_count[PYRAMID]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Prisms", &type_count[PRISM]);
    HDF_WriteScalar(h5out[i], directoryBase, "Number Of Hexahedron", &type_count[HEX]);
  }
  
  ////////////////////////////////////////////////////////
  //
  //  NODE DATA WRITING
  //
  ////////////////////////////////////////////////////////

  //write ghost nodes owning process ids
  for(i = 0; i < np; i++){
    for(j = 0; j < ghostNodesCount[i]; j++){
      Int node = ghostNodes[i*goffset + j];
      ghostDoms[j] = part[node];
    }
    directoryBase = "/Mesh/";
    HDF_WriteArray(h5out[i], directoryBase, "Ghost Nodes Owning Process", 
		   &ghostDoms[0], ghostNodesCount[i]);
  }
  //build hash_map to lookup local id of nodes on owning processes
  //also used to reference local elems by local node number system
  Int* hash_nodes = new Int[nnode];
  for(i = 0; i < np; i++){
    for(j = 0; j < localNodesCount[i]; j++){
      Int node = localNodes[localNodesOffset[i] + j];
      //set node index to location in whatever process'
      //local list it is currently in 
      hash_nodes[node] = j;
    }
  }
  //write local id of ghost nodes
  for(i = 0; i < np; i++){
    vint.clear();
    directoryBase = "/Mesh/";
    for(j = 0; j < ghostNodesCount[i]; j++){
      Int node = ghostNodes[i*goffset + j];
      vint.push_back(hash_nodes[node]);
    }
    if(vint.size() != (UInt)ghostNodesCount[i]){
      std::cerr << "WARNING vint size does not match ghost nodes count" << std::endl;
    }
    HDF_WriteArray(h5out[i], directoryBase, "Ghost Nodes Local Id", 
		   &vint.front(), vint.size());
  }
  //write local and ghost nodes coords
  for(i = 0; i < np; i++){
    vreal.clear();
    directoryBase = "/Mesh/";
    for(j = 0; j < localNodesCount[i]; j++){
      Int node = localNodes[localNodesOffset[i] + j];
      vreal.push_back(xyz[node*3 + 0]);
      vreal.push_back(xyz[node*3 + 1]);
      vreal.push_back(xyz[node*3 + 2]);
    }
    for(j = 0; j < ghostNodesCount[i]; j++){
      Int node = ghostNodes[i*goffset + j];
      vreal.push_back(xyz[node*3 + 0]);
      vreal.push_back(xyz[node*3 + 1]);
      vreal.push_back(xyz[node*3 + 2]);
    }
    HDF_WriteArray(h5out[i], directoryBase, "Nodal Coordinates", 
		   &vreal.front(), vreal.size());
  }

  ////////////////////////////////////////////////////////
  //
  //  ELEMENT DATA WRITING
  //
  ////////////////////////////////////////////////////////
  
  //write local and split elems factags
  for(i = 0; i < np; i++){
    vint.clear();
    directoryBase = "/Mesh/";
    for(j = 0; j < localElemsCount[i]; j++){
      Int gelem = localElems[i][j];
      Element<Real>& element = *m.elementList[gelem];
      Int factag = element.GetFactag();
      vint.push_back(factag);
    }
    for(j = 0; j < splitElemsCount[i]; j++){
      Int gelem = splitElems[i*seoffset + j];
      Element<Real>& element = *m.elementList[gelem];
      Int factag = element.GetFactag();
      vint.push_back(factag);
    }
    HDF_WriteArray(h5out[i], directoryBase, "Element Factags", 
		   &vint.front(), vint.size());
  }
  //write local and split elem data (local numbering)
  for(i = 0; i < np; i++){
    vint.clear();
    directoryBase = "/Mesh/";
    for(j = 0; j < localElemsCount[i]; j++){
      Int gelem = localElems[i][j];
      Element<Real>& element = *m.elementList[gelem];
      Int etype = element.GetType();
      Int* nodes; 
      Int nnodes = element.GetNodes(&nodes);
      vint.push_back(etype);
      for(k = 0; k < nnodes; k++){
	Int node = nodes[k];
      	//translate node from global number to number local to current process
	//node is guaranteed to be local from the way we defined local elems
	node = hash_nodes[node];
	vint.push_back(node);
      }
    }
    for(j = 0; j < splitElemsCount[i]; j++){
      Int gelem = splitElems[i*seoffset + j];
      Element<Real>& element = *m.elementList[gelem];
      Int etype = element.GetType();
      Int* nodes; 
      Int nnodes = element.GetNodes(&nodes);
      vint.push_back(etype);
      for(k = 0; k < nnodes; k++){
	Int node = nodes[k];
	Int* node2;
	Int pos;
	//translate node from global number to number local to current process
	//node is not guaranteed to be in local list but will be in the union
	//of the local list, and ghost list
	if(part[node] == i){
	  //if node is local read through hash map
	  node = hash_nodes[node];
	}
	else{
	  //if node is not local it must be in the ghost list
	  //if this fails something is really wrong above
	  node2 = find(&ghostNodes[i*goffset], &ghostNodes[i*goffset] + ghostNodesCount[i], node);
	  if(node2 == &ghostNodes[i*goffset] + ghostNodesCount[i]){
	    cerr << "Cannot translate node " << node << " from global to local on process " << i << endl;
	  }
	  pos = node2 - &ghostNodes[i*goffset];
	  node = localNodesCount[i] + pos;
	}
	vint.push_back(node);
      }
    }
    HDF_WriteArray(h5out[i], directoryBase, "Element Data", 
		   &vint.front(), vint.size());
  }

  for(i = 0; i < np ; i++){
    HDF_CloseFile(h5out[i]);
  }

  delete [] hash_nodes;
  delete [] localNodes;
  delete [] ghostNodes;
  delete [] ghostDoms;
  for(i = 0; i < np; i++){
    delete [] localElems[i];
  }
  delete [] localElems;
  delete [] splitElems;
  
  delete [] part;
  
  return (ierr);
}

void Coordpartition(Mesh<Real>* m, Int* part, Int np, Int* nodesPart)
{
  //split into np divisions based on x-coords
  Int nnode = m->GetNumNodes();
  double const * xyz = m->GetNodeCoords();
  double maxX = xyz[0*3 + 0];
  double minX = xyz[0*3 + 0];
  for(Int i = 1; i < nnode; i++){
    minX = MIN(minX, xyz[i*3 + 0]);
    maxX = MAX(maxX, xyz[i*3 + 0]);
  }
  double span = maxX - minX;
  double div[np+1];
  div[0] = minX - 0.01;
  cout << "DIVISIONS: " << endl;
  cout << "=======================" << endl;
  for(Int i = 1; i < np; i++){
    div[i] = div[i-1] + span/(double)np + 0.1;
  }
  //to avoid roundoff missing points
  div[np] = maxX + 0.01;
  for(Int i = 0; i < np; i++){
    cout << "Partition " << i << ": " << div[i] << " to " << div[i+1] << endl;
  }
  cout << endl;

  for(Int j = 0; j < np; j++){
    double low = div[j];
    double high = div[j+1];
    for(Int i = 0; i < nnode; i++){
      double x = xyz[i*3 + 0];
      if(x > low && x < high){
	part[i] = j;
      }
    }
  }
  for(Int i = 0; i < nnode; i++){
    nodesPart[part[i]]++;
  }
  cout << "COORDX: Nodes per partition" << endl;
  cout << "==========================" << endl;
  for(Int i = 0; i < np; i++){
    cout << "\tPartition " << i << ": " << nodesPart[i] << endl;
  }

}

void METISpartition(Mesh<Real>* m, Int* part, Int np, Int* nodesPart)
{
  Int i;
  
  //calculate weights for vertics based on number of connections
  Int nnode = m->GetNumNodes();
  Int* edgesw = NULL;
  Int* vertw = new Int[nnode];
  Int* vertsize = new Int[nnode];
  Int edgecut;
  real_t* ubvec = NULL;
  real_t* tpwgts = NULL;
  idx_t ncon = 1;
  //use default options
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  for(i = 0; i < nnode; i++){
    Int order = m->ipsp[i+1] - m->ipsp[i];
    vertw[i] = order;
    vertsize[i] = order;
  }
  Int partWeight[np];
  Int sum;
  for(i = 0; i < np; i++){
    nodesPart[i] = 0;
    partWeight[i] = 0;
  }
  
  if(np == 1){
    cout << "METIS: NP = 1 -- No partitioning required" << endl;
    for(i = 0; i < nnode; i++){
      part[i] = 0;
    }
    edgecut = 0; 
  }
  
  else{
    cout << "METIS: Partitioning graph" << endl;
    if(np < 8){
      Int lnnode = nnode;
      METIS_PartGraphRecursive(&lnnode, &ncon, (idx_t*)m->ipsp, (idx_t*)m->psp, (idx_t*)vertw, (idx_t*)vertsize,
			       (idx_t*)edgesw, &np, tpwgts, ubvec, options, &edgecut, (idx_t*)part); 
    }
    else{
      Int lnnode = nnode;
      METIS_PartGraphKway(&lnnode, &ncon, (idx_t*)m->ipsp, (idx_t*)m->psp, (idx_t*)vertw, (idx_t*)vertsize, 
			  (idx_t*)edgesw, &np, tpwgts, ubvec, options, &edgecut, (idx_t*)part); 
      
    }
  }
  cout << "METIS: Edges cut: " << edgecut << endl;
  for(i = 0; i < nnode; i++){
    nodesPart[part[i]]++;
  }
  cout << endl;
  cout << "METIS: Nodes per partition" << endl;
  cout << "==========================" << endl;
  for(i = 0; i < np; i++){
    cout << "\tPartition " << i << ": " << nodesPart[i] << endl;
  }
  cout << endl;
  cout << "METIS: Unweighted load balancing" << endl;
  cout << "================================" << endl;

  for(i = 0; i < np; i++){
    cout << "\tPartition " << i << ": " << (Real)nodesPart[i]/(Real)nnode << endl;
  }
  cout << endl;
  cout << "METIS: Weighted load balancing" << endl;
  cout << "==============================" << endl;
  sum = 0;
  for(i = 0; i < nnode; i++){
    Int weight = m->ipsp[i+1] - m->ipsp[i];
    sum += weight;
    partWeight[part[i]] += weight;
  }
  for(i = 0; i < np; i++){
    cout << "\tPartition " << i << ": " << (Real)partWeight[i]/(Real)sum << endl;
  }
  cout << endl;

  
}
