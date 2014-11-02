#include "io.h"

namespace STRUCTDYN{

int WriteVTK_Ascii(std::string casename, int* nelem, Element* elems, int nnode, 
		   double* xyz, double* variables, int nvars, std::string* varNames)
{
  int err = 0;

  std::ofstream fout;
  int i, j, k, e;
  std::string filename = casename + ".vtk";
  int elemtypes = 1;
  int totalElems = 0;
  int count = 0;

  std::cout << "VTK ASCII I/O: Writing solution file --> " << filename << std::endl;

  fout.open(filename.c_str());

  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << "Solver solution data" << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  fout << "POINTS " << nnode << " double" << std::endl;
  for(i = 0; i < nnode; i++){
    fout << xyz[i*3 + 0] << " " << xyz[i*3 + 1] << " " << xyz[i*3 + 2] << std::endl;
  }

  for(i = 0; i < elemtypes; i++){
    totalElems += nelem[i];
  }

  for(i = 0; i < totalElems; i++){
    count += elems[i].GetNnodes();
  }

  fout << "CELLS " << totalElems << " " << totalElems+count << std::endl;
  for(i = 0; i < totalElems; i++){
    fout << elems[i].GetNnodes() << " "; 
    for(j = 0; j < elems[i].GetNnodes(); j++){ 
      fout << elems[i].nodes[j] << " ";
    }
    fout << "\n";
  }

  fout << "CELL_TYPES " << totalElems << std::endl;
  for(i = 0; i < elemtypes; i++){
    for(j = 0; j < nelem[i]; j++){
      switch (i) {
	//1D beams
      case 0 :
	fout << "3" << std::endl;
	break;
      default:
	std::cerr << "Type not defined WriteVTK_Ascii()" << std::endl;
      }
    }
  }

  //if we have cell data for some reason do this... 
  //fout << "CELL_DATA " << totalElems << std::endl;
  //fout << "SCALARS BoundaryIds int " << 1 << std::endl;

  fout << "POINT_DATA " << nnode << std::endl;
  for(i = 0; i < nvars; i++){
    fout << "SCALARS " << varNames[i] << " double " << 1 << std::endl;
    fout << "LOOKUP_TABLE default" << std::endl;
    for(j = 0; j < nnode; j++){
      fout << variables[i*nnode + j] << "\n";
    }
  }
  
  fout << std::endl;
  fout.close();
  std::cout << "VTK ASCII I/O: File write successful!!" << std::endl;


  return err;
}

}
