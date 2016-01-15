#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "mesh.h"
#include "general.h"


using namespace std;

int main(int argc, char* argv[]){

  Int i;
  Int pos;
  Int ierr = 0;
  Mesh<Real> m;
  string filename;
  string casename;
  string fileextension;
  stringstream iss;
  Real txyz[3];
  Real tol;
  Int ptid;

  if(!((argc == 6) || (argc == 3))){
    cerr << "Invalid arguments!!!" << endl;
    cerr << "USE: " << argv[0] <<" <gridname> <xpos> <ypos> <zpos> <tolerance>" << endl;
    return (1);
  }

  filename = argv[1];
  if(argc == 6){
    iss.clear();
    iss.str(argv[2]);
    iss >> txyz[0];
    iss.clear();
    iss.str(argv[3]);
    iss >> txyz[1];
    iss.clear();
    iss.str(argv[4]);
    iss >> txyz[2];
    iss.clear();
    iss.str(argv[5]);
    iss >> tol;
  }
  else if(argc == 3){
    iss.clear();
    iss.str(argv[2]);
    iss >> ptid;
  }

  pos = filename.rfind(".");
  casename = filename.substr(0,pos);
  fileextension = filename.substr(pos+1, filename.size() - (casename.size() ));
  if(fileextension == "crunch"){
    m.ReadCRUNCH_Ascii(filename);
  }
  else if(fileextension == "ugrid"){
    m.ReadUGRID_Ascii(filename);
  }
  else{
    cerr << "File extension " << fileextension << " not recognized" << endl;
    return(2);
  }

  if(argc == 6){
    cout << "Scanning for points within " << tol << " of " << txyz[0] 
	 << " " << txyz[1] << " " << txyz[2] << endl << endl;
    //scan x, then y, then z
    for(i = 0; i < m.nnode; i++){
      if(fabs(m.xyz[i*3 + 0] - txyz[0]) < tol){
	if((fabs(m.xyz[i*3 + 1] - txyz[1]) < tol) && (fabs(m.xyz[i*3 + 2] - txyz[2]) < tol) ){
	  cout << "ptid: " << i << " within tolerance -- " << m.xyz[i*3 + 0] << " " 
	       << m.xyz[i*3 + 1] << " " << m.xyz[i*3 + 2] <<  endl;
	}
      }
    }
  }
  else if(argc == 3){
    cout << "Printing coordinates for point: " << ptid << endl;
    cout << m.xyz[ptid*3 + 0] << " " << m.xyz[ptid*3 + 1] << " " << m.xyz[ptid*3 + 2] <<  endl;
	
  }

  return(ierr);
}
