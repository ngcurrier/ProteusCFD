#include "mesh.h"
#include "general.h"
#include "h5layer.h"
#include "mem_util.h"
#include "macros.h"
#include "boundaryInflation.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  HDF_TurnOffErrorHandling();
  MPI_Init(&argc, &argv);
  PObj<Real> pobj;

  //mesh has not been reordered at this point
  Bool reordered = false;

  Mesh<Real> m;
  m.SetParallelPointer(&pobj);
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

  if(argc != 2){
    cerr << "Invalid arguments" << endl;
    cerr << argv[0] << " <Mesh file>" << endl;
    return(3);
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
    ierr = m.ReadGMSH_Master(filename);
  }
  else if(fileextension == "cgns"){
#ifdef _HAS_CGNS
    ierr = m.ReadCGNS(filename);
#else
    Abort << "Proteus was not compiled with CGNS support. Please rebuild.";
#endif
  }
  else{
    cerr << "File extension " << fileextension << " not recognized" << endl;
    return 0;
  }

  if(ierr){
    cerr << "Grid read failure on file -- " << filename << endl;
    return 1;
  }

  //grab node coords
  double const * xyz = m.GetNodeCoords();
  
  //build utility maps for METIS
  m.BuildMaps();
  m.CalcMetrics();
  pobj.BuildCommMaps(&m);

  Int nnode = m.GetNumNodes();
  Int lnelem = m.GetNumLocalElem();

  std::string tags;
  std::cout << "Please list the boundary ids to generate boundary layers on separated by commas:" << std::endl;
  std::getline(std::cin, tags);
  
  std::vector<std::string> boundaries = Tokenize(tags, ',');
  std::vector<Real> boundaryThicknesses;
  std::vector<int> numberOfLayers;
  std::vector<int> boundaryFactagList;

  for(int i = 0; i < boundaries.size(); ++i){
    std::stringstream ss(boundaries[i]);
    int tag;
    ss >> tag;
    boundaryFactagList.push_back(tag);
    std::cout << "For tag " << tag << " please input layer thickness: " << std::endl;
    std::string tmp;
    std::getline(std::cin, tmp);
    ss.clear();
    ss.str(tmp);
    Real thickness;
    ss >> thickness;
    boundaryThicknesses.push_back(thickness);
    std::cout << "For tag " << tag << " please input number of cell layers: " << std::endl;
    std::getline(std::cin, tmp);
    ss.clear();
    ss.str(tmp);
    int nlayers;
    ss >> nlayers;
    numberOfLayers.push_back(nlayers);
  }

  std::cout << "Boundary inflation parameters:" << std::endl;
  for(int i = 0; i < boundaryFactagList.size(); ++i){
    std::cout << boundaryFactagList[i] << ":\t" << boundaryThicknesses[i] << "\t" << numberOfLayers[i] << std::endl;
  }

  std::cout << "Generating boundary layers" << std::endl;
  Real growthRate = 1.2;
  GenerateBoundaryLayers(boundaryFactagList, boundaryThicknesses, numberOfLayers, &m, growthRate);

  //write out the modified grid
  std::string newGridFilename = casename + ".new";
  m.WriteCRUNCH_Ascii(newGridFilename);

  MPI_Finalize();
}
