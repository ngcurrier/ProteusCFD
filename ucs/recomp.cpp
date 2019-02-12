#include "mesh.h"
#include "etypes.h"
#include "general.h"
#include "strings_util.h"
#include "h5layer.h"
#include "solutionField.h"
#include "boundaryInflation.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm> //std::sort
#include <set>

using namespace std;


int main(int argc, char* argv[])
{
  Int np = 0;
  stringstream temposs;
  string casename;
  string temps;
  Int format = 0;
  Mesh<Real> m;
  
  if(argc != 2 && argc != 3 ){
    cerr << "Invalid arguments" << endl;
    cerr << argv[0] << " <casename> <format flag>" << endl;
    return (1);
  }
  
  casename = argv[1];
  if(argc == 3){
    temposs.clear();
    temposs << argv[2];
    temposs >> format;
    temposs.clear();
  }
  HDF_TurnOffErrorHandling();
  hid_t h5temp = HDF_OpenFile(casename + ".0.h5", 0);
  if(h5temp < 0){
    std::cerr << "RECOMP I/O: Could not open file -- " << casename + ".0.h5" << std::endl;
    return -1;
  }
  //
  //  read standard data from HDF files across all processes
  //  
  
  Int reordered, rescaled;
  HDF_ReadScalar(h5temp, "/", "Number Of Processors", &np);
  //read flag to tell if mesh has been previously reordered
  HDF_ReadScalar(h5temp, "/", "Reordered", &reordered);
  if(reordered){
    m.SetMeshReordered();
  }
  HDF_ReadScalar(h5temp, "/", "Rescaled", &rescaled);
  if(rescaled){
    m.SetMeshScaled();
  }
  //read total number of nodes in mesh
  HDF_ReadScalar(h5temp, "/", "Global Number Of Nodes", &m.gnnode);
  //read total number of elems in mesh
  HDF_ReadScalar(h5temp, "/", "Global Number Of Elements", &m.gnelem);
  //read number of factags in mesh
  HDF_ReadScalar(h5temp, "/", "Number Of Factags", &m.nfactags);

  //we will read in the global mesh, set local values appropriately
  m.lnelem = m.gnelem;
  m.nnode = m.gnnode;
  m.gnode = 0;
  m.MemInitMesh();

  HDF_CloseFile(h5temp);

  //now open the file series necessary to obtain all the data
  hid_t h5in[np];
  cout << "RECOMP I/O: Reading " << np << " partitioned files" << endl;
  cout << "RECOMP I/O: Opening files -- " << casename + ".*.h5" << endl;
  for(Int i = 0; i < np; i++){
    temposs.str("");
    temposs << i;
    string filename = casename + "." + temposs.str() + ".h5";
    cout << "RECOMP I/O: Opening file -- " << filename << endl;
    h5in[i] = HDF_OpenFile(filename, 0);
    if(h5in[i] < 0){
      std::cerr << "RECOMP I/O: Could not open file -- " << filename << std::endl;
      return -1;
    }
  }

  //TODO: test for unsteady grid movement, needed for moving writes

  Int nnode[np];
  Int gnode[np];
  Int localElemCount[np];
  Int nodeOffset[np];
  for(Int i = 0; i < np; i++){
    Int temp = 0;
    localElemCount[i] = 0;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Local Nodes", &nnode[i]);
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Ghost Nodes", &gnode[i]);
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Triangles", &temp);
    localElemCount[i] += temp;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Quadrilaterals", &temp);
    localElemCount[i] += temp;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Tetrahedron", &temp);
    localElemCount[i] += temp;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Pyramids", &temp);
    localElemCount[i] += temp;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Prisms", &temp);
    localElemCount[i] += temp;
    HDF_ReadScalar(h5in[i], "/Mesh/", "Number Of Hexahedron", &temp);
    localElemCount[i] += temp;
  }
  nodeOffset[0] = 0;
  for(Int i = 1; i < np; i++){
    nodeOffset[i] = nodeOffset[i-1] + nnode[i-1];
  }

  //read in xyz coordinates
  for(Int i = 0; i < np; i++){
    cout << "RECOMP I/O: Reading coordinates from process " << i << endl;
    Int ncoords = (nnode[i]+gnode[i])*3;
    Real* xyztemp = new Real[ncoords];
    HDF_ReadArray(h5in[i], "/Mesh/", "Nodal Coordinates", &xyztemp, &ncoords);
    memcpy(&m.xyz[nodeOffset[i]*3], xyztemp, sizeof(Real)*nnode[i]*3);
    delete [] xyztemp;
  }

  //read in element list on each partitioned mesh
  std::vector<Element<Real>*> elementList;
  //use to get a unique sorted list later
  std::set<Element<Real>*, ElementPtrComp<Real> > splitElementList;
  Int totalElemEstimate = 0;
  for(Int i = 0; i < np; i++){
    totalElemEstimate += localElemCount[i];
  }
  elementList.reserve(totalElemEstimate);
  for(Int i = 0; i < np; i++){
    cout << "RECOMP I/O: Reading elements from process " << i << endl;
    Int* factagTempLocal = new Int[localElemCount[i]];
    Int* gnodeOwningProcess = new Int[gnode[i]];
    Int* gnodeLocalId = new Int[gnode[i]];
    Int listSize = -1;
    Int* elementData = NULL;
    //get size of array
    HDF_ReadArray(h5in[i], "/Mesh/", "Element Data", &elementData, &listSize);
    elementData = new Int[listSize];
    //read arrays
    HDF_ReadArray(h5in[i], "/Mesh/", "Element Data", &elementData, &listSize);
    HDF_ReadArray(h5in[i], "/Mesh/", "Element Factags", &factagTempLocal, &localElemCount[i]);	  
    HDF_ReadArray(h5in[i], "/Mesh/", "Ghost Nodes Local Id", &gnodeLocalId, &gnode[i]);
    HDF_ReadArray(h5in[i], "/Mesh/", "Ghost Nodes Owning Process", &gnodeOwningProcess, &gnode[i]);
    //read in element data and put factags where they belong
    //they are stored as [type-e1, node1-e1, node2-e1, .... type-e2] in a linear array
    Int loc = 0;
    for(Int j = 0; j < localElemCount[i]; j++){
      Bool isSplit = false;
      Int nodes[8]; 
      Int e = elementData[loc];
      loc++;
      //increment the node id by the nodeoffset value
      for(Int k = 0; k < m.GetNumElemNodes(e); k++){
	Int localNodeId = elementData[loc];
	//check if the element node is a ghost node
	//if so, we have to look up that ghost node's global id
	if(localNodeId >= nnode[i]){
	  isSplit = true;
	  Int gnodeId = localNodeId - nnode[i];
	  Int owner = gnodeOwningProcess[gnodeId];
	  Int localId = gnodeLocalId[gnodeId];
	  nodes[k] = localId + nodeOffset[owner];
	}
	//node is local use standard offset
	else{
	  nodes[k] = localNodeId + nodeOffset[i];
	}
	loc++;
      }
      Element<Real>* tempe = NULL;
      switch (e) {  
      case TRI :
	tempe = new Triangle<Real>;
	break;
      case QUAD :
	tempe = new Quadrilateral<Real>;
	break;
      case TET :
	tempe = new Tetrahedron<Real>;
	break;
      case PYRAMID :
	tempe = new Pyramid<Real>;
	break;
      case PRISM :
	tempe = new Prism<Real>;
	break;
      case HEX :
	tempe = new Hexahedron<Real>;
	break;
      default :
	std::cerr << "Type not defined in ReadPartedMesh()" << std::endl;
      }
      tempe->Init(nodes);
      tempe->SetFactag(factagTempLocal[j]);
      Int oldsize = splitElementList.size();
      if(isSplit){
	splitElementList.insert(tempe);
	if(splitElementList.size() != oldsize){
	  m.nelem[e]++;
	}
      }
      else{
	elementList.push_back(tempe);
	m.nelem[e]++;
      }
    }
    delete [] elementData;
    delete [] factagTempLocal;
    delete [] gnodeOwningProcess;
    delete [] gnodeLocalId;
  }

  //sort the list and then delete duplicates
  m.elementList.reserve(splitElementList.size()+elementList.size());
  for(std::vector<Element<Real>*>::iterator it = elementList.begin(); it != elementList.end(); ++it){
    Element<Real>* element = *it;
    m.elementList.push_back(element);
  }
  for(std::set<Element<Real>*>::iterator it = splitElementList.begin(); 
      it != splitElementList.end(); ++it){
    m.elementList.push_back(*it);
  }
  cout << "RECOMP I/O: Read " << splitElementList.size()+elementList.size() 
       << " elements from file" << endl;
  //clear lists, pointers now stored by the mesh
  elementList.clear();
  splitElementList.clear();

  if(m.lnelem != m.elementList.size()){
    std::cerr << "WARNING: expected number of elements is " << m.lnelem << " found " 
	      << m.elementList.size() << std::endl;
    m.lnelem = m.elementList.size();
    //return (-1);
  }

  std::vector<int> boundaryFactagList;
  std::vector<Real> boundaryThicknesses;
  std::vector<int> numberOfLayers;
  GenerateBoundaryLayers(boundaryFactagList, boundaryThicknesses, numberOfLayers, &m);
  
  std::vector<SolutionField<Real>*> fields;

  //Read in the datasets from one of the files
  std::vector<std::string> dataPaths;
  HDF_GetAllDatasetsInFile(h5in[0], dataPaths);

  //check to see if the datasets are in a solution directory
  Bool solutionPresent = false;
  std::vector<std::string> solutionVect;

  for(std::vector<std::string>::iterator it = dataPaths.begin(); it != dataPaths.end(); ++it){
    std::string& path = *it;
    std::string keyword = "/Solution/";
    size_t loc = path.find(keyword);
    if(loc != std::string::npos){
      solutionPresent = true;
      solutionVect.push_back(*it);
    }
  }
  //create a list which identifies the solution set of each entry, this lets us write whole solution sets
  //instead of just each entry
  Int sets[solutionVect.size()+1];
  UInt count = 0;
  sets[0] = count;
  Int s = 0;
  std::string first;
  if(solutionVect.size() != 0){
    first = *solutionVect.begin();
  }
  std::string oldvalue = first.substr(0, first.rfind("/")+1);
  for(std::vector<std::string>::iterator it = solutionVect.begin(); it != solutionVect.end(); ++it){
    std::string& path = *it;
    size_t loc = path.rfind("/");
    std::string newvalue = path.substr(0, loc+1);
    if(newvalue != oldvalue){
      count++;
      oldvalue = newvalue;
    }
    sets[s] = count;
    s++;
  }
  sets[solutionVect.size()] = count;
  std::string solfilename = casename;
  if(solutionPresent == false){
    std::cout << "RECOMP I/O: Solution not present, have you run the solver yet? ;) " << std::endl;
    if(format == 0){
      m.WriteXMLVTK_Binary(solfilename, fields);
    }  
    else if(format == 1){
      m.WriteVTK_Binary(solfilename, fields);
    }
    else if(format == 2){
      m.WriteCRUNCH_Ascii(solfilename);
    }
    else if(format == 3){
      m.WriteVTK_Ascii(solfilename, fields);
    }
    else{}
  }
  else{ //if solution is present check all the fields for writing, etc.

    // if we are writing XML based vtk the also create a pvd file for including
    // time stepping context in the visualization
    std::ofstream fpvd;
    if(format == 0){
      std::string sscase;
      sscase = casename + ".pvd";
      fpvd.open(sscase.c_str());
      fpvd << "<?xml version=\"1.0\"?>" << std::endl;
      fpvd << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
      fpvd << "<Collection>";
    }
    
    //loop over all solutions that are present (i.e. a time series), recompose and write one at a time
    count = 0;
    for(std::vector<std::string>::iterator it = solutionVect.begin(); it != solutionVect.end(); ++it){
      double time = 0.0;
      std::string& path = *it;
      std::string timestepstring = "timestep-";
      std::string stepnumber = "";
      solfilename = casename;
      size_t loc = path.rfind('/')+1;
      std::string directory = path.substr(0, loc);
      std::string dataname = path.substr(loc, std::string::npos);
      size_t loct = path.find(timestepstring);
      //if timestepping is in file, parse the step level
      if(loct != std::string::npos){
	//find the next separator
	size_t locsep = path.find('/', loct);
	if(locsep == std::string::npos){
	  std::cerr << "WARNING: Timestep path malformed -- " << path << std::endl;
	}
	size_t length = timestepstring.size();
	loct += length;
	stepnumber = path.substr(loct,locsep-loct);
	//add the timestep number to the file base name
	solfilename += "-" + stepnumber;
	HDF_ReadScalar(h5in[0], "/SolutionTime/timestep-" + stepnumber + "/", "time", &time);
      }
      SolutionField<Real>* field = new SolutionField<Real>(m, h5in[0], path);
      fields.push_back(field);
      //the field exists and memory is allocated, now read in the parallel partitioned data
      //and reassemble it in a meaningful way
      for(Int i = 0; i < np; i++){
	Int size = -1;
	Real* tempdata = NULL;
	//probe dataset for size
	HDF_ReadArray(h5in[i], directory, dataname, &tempdata, &size);
	tempdata = new Real[size];
	//read the dataset
	HDF_ReadArray(h5in[i], directory, dataname, &tempdata, &size);
	//copy the appropriate chunk to the field, i.e. the first nnode values
	Real* fielddata = field->GetData(FIELDS::STATE_NONE);
	Int ndof = field->GetNdof();
	memcpy(&fielddata[nodeOffset[i]*ndof], tempdata, sizeof(Real)*nnode[i]*ndof);
	delete [] tempdata;
      }
      //check for time sensitive coords., if present replace current values for moving mesh 
      std::string xyzDirectory = "/Mesh/" + timestepstring + stepnumber +"/";
      if(HDF_TestDirectory(h5in[0], xyzDirectory)){
	//read in xyz coordinates
	for(Int i = 0; i < np; i++){
	  Int ncoords = (nnode[i]+gnode[i])*3;
	  Real* xyztemp = new Real[ncoords];
	  HDF_ReadArray(h5in[i], xyzDirectory, "Nodal Coordinates", &xyztemp, &ncoords);
	  memcpy(&m.xyz[nodeOffset[i]*3], xyztemp, sizeof(Real)*nnode[i]*3);
	  delete [] xyztemp;
	}
      }

      if((sets[count] != sets[count+1]) || count == solutionVect.size()-1){
	if(format == 0){
	  m.WriteXMLVTK_Binary(solfilename, fields);
	  fpvd << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << solfilename+".vtu" << "\"/>" << std::endl;
	}  
	else if(format == 1){
	  m.WriteVTK_Binary(solfilename, fields);
	}
	else if(format == 2){
	  m.WriteCRUNCH_Ascii(solfilename);
	}
	else if(format == 3){
	  m.WriteVTK_Ascii(solfilename, fields);
	}
	else{}
	//clear the fields for the next pass
	for(std::vector<SolutionField<Real>*>::iterator it = fields.begin(); it != fields.end(); ++it){
	  delete *it;
	}
	fields.clear();
      }
      count++;
    }
    if(format == 0){
      fpvd << "</Collection>" << std::endl;
      fpvd << "</VTKFile>" << std::endl;
      fpvd.close();
    }
  } // end solution writing loop

  
  //close hdf file series
  for(Int i = 0; i < np; i++){
    HDF_CloseFile(h5in[i]);
  }

  return 0;
}
