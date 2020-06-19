#ifdef _HAS_PYTHON
// we have to include this first b/c python doesn't check for _XOPEN_SOURCE
// and hasn't fixed it yet
#include <Python.h>
#endif

#include "mesh.h"
#include "general.h"
#include "h5layer.h"
#include "mem_util.h"
#include "macros.h"
#include "boundaryInflation.h"
#include "bc.h"
#include <unistd.h>
#include <getopt.h>
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
  BoundaryConditions<Real>* bc = NULL;

  //mesh has not been reordered at this point
  Bool reordered = false;

  Mesh<Real> m;
  m.SetParallelPointer(&pobj);
  string filename;
  string outfilename;
  string fileextension;
  string casename;
  Real scale = 1.0;
  size_t pos;
  Int np = 0;
  
  //use these as dynamic holders, fast and easy
  vector<Int> vint;
  vector<Real> vreal;

  Int i, j, k;
  Int ierr = 0;

  
  int opt; 
  int sfnd = false;
  bool insertBoundaryLayer = false;
  
  // put ':' in the starting of the 
  // string so that program can  
  //distinguish between '?' and ':'  
  while(1){
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"input",          required_argument, 0,  'i'},
      {"help",           no_argument,       0,  'h'},
      {"scale",          required_argument, 0,  's'},
      {"boundarylayer",  no_argument,       0,  'b'},
      {"output",         required_argument, 0,  'o'},
      {0,         0,                 0,  0 }
    };

    // one colon following option means the option takes an arg
    // two colons following option means the option takes an optional arg
    opt = getopt_long(argc, argv, "s:i:o:hb",
		    long_options, &option_index);
    if (opt == -1)
      break;

    std::string tmp;
    std::stringstream ss;
    switch(opt)  
      {  
      case 'i':
	if (!optarg){
	  cerr << "WARNING: input file option requires a filename" << std::endl;
	  return 3;
	}
	cout << "Input filename: " << optarg << std::endl;
	filename = optarg;
	break;
      case 'o':
	if (!optarg){
	  cerr << "WARNING: output file option requires a filename" << std::endl;
	  return 3;
	}
	cout << "Output filename: " << optarg << std::endl;
	outfilename = optarg;
	break;
      case 's':
	if (!optarg){
	  cerr << "WARNING: require scale factor" << std::endl;
	  return 3;
	}
	tmp = optarg;
	ss.clear();
	ss << tmp;
	ss >> scale;
	cout << "Scale factor: " << scale << std::endl;
	sfnd = true;
	break;
      case 'b':
	insertBoundaryLayer = true;
	break;
      case ':':  
	printf("option needs a value\n");  
	break;
      case 'h':
	cout << "USAGE: " << argv[0] << " <options>" << std::endl;
	cout << "Options:" << std::endl;
	cout << "--------------------------------------------------------" << std::endl;
	cout << "--input\t\t-i\t<Input file name>" << std::endl;
	cout << "--output\t-o\t<Output file name>" << std::endl;
	cout << "--boundarylayer\t-b" << std::endl;
	cout << "--scale\t\t-s\t<Scale factor for mesh, i.e. reference length desired>" << std::endl;
	return 0;
	break;
      default:
	cout << "WARNING: unknown option: " << static_cast<char>(optopt) << std::endl;
	return 3;
	break;  
      }  
    }  
  
  // optind is for the extra arguments 
  // which are not parsed 
  for(; optind < argc; optind++){
    cout << "WARNING: option not parsed \'" << argv[optind] << "\'" << std::endl;
  }
  

  if (filename.size() == 0){
    cerr << "WARNING: filename option -i not captured" << endl;
    return 3;
  }
  
  //input file name parsing -- reader selection
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

  if(insertBoundaryLayer){
    // Read boundary condition file casename.bc
    bc = new BoundaryConditions<Real>;
    Int ierr = bc->ReadFile(casename, m.GetMaximumFactag());
    if(ierr){
      Abort << "BC file read failed";
    }

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
      std::cout << "For tag " << tag << " please input first layer thickness: " << std::endl;
      std::string tmp;
      std::getline(std::cin, tmp);
      ss.clear();
      ss.str(tmp);
      Real thickness;
      ss >> thickness;
      boundaryThicknesses.push_back(thickness);
      std::cout << "For tag " << tag << " please input number of cell layers <negative value matches grid spacing>: " << std::endl;
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
    // this is the maximum ecommended value by the DPW gridding guidelines
    Real growthRate = 1.25;
    GenerateBoundaryLayers(boundaryFactagList, boundaryThicknesses, numberOfLayers, bc, &m, growthRate);
  }
  
  //this does a direct scaling of the mesh by this value (i.e. direct multiplication)
  if(sfnd){
    m.ScaleMesh(1.0/scale);
  }
  
  //write out the modified grid
  pos = outfilename.rfind(".");
  casename = outfilename.substr(0,pos);
  fileextension = outfilename.substr(pos+1, outfilename.size() - (casename.size() ));
  if(fileextension == "crunch"){
    m.WriteCRUNCH_Ascii(casename);
  }
  else if(fileextension == "ugrid"){
    m.WriteUGRID_Ascii(casename);
  }
  else if(fileextension == "vtk"){
    std::vector<SolutionField<Real>*> fields;
    m.WriteVTK_Ascii(casename, fields);
  }
  else if(fileextension == "stl"){
    std::cout << "\nInput factag to extract (negative extracts all): " << std::endl;
    std::stringstream ss2;
    std::string tmp2;
    std::getline(std::cin, tmp2);
    ss2.str(tmp2);
    Int factag;
    ss2 >> factag;
    m.WriteSTL_Ascii(casename, factag);
  }
  else if(fileextension == "cas"){
    // Read boundary condition file casename.bc
    bc = new BoundaryConditions<Real>;
    Int ierr = bc->ReadFile(casename, m.GetMaximumFactag());
    if(ierr){
      std::cout << "BC file read failed, will not port BCs to FLUENT case file";
    }
    m.WriteFluentCase_Ascii(casename, bc);
  }
  else{
    cerr << "File extension " << fileextension << " not recognized for output mesh" << endl;
    return 0;
  }


  MPI_Finalize();
}
