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
#include <vector>

using namespace std;

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
  m.BuildMapsDecomp();

  Int nnode = m.GetNumNodes();
  Int lnelem = m.GetNumLocalElem();


}
