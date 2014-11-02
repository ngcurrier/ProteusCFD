#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "portFunctions.h"
#include "param.h"
#include "portGlobal.h"
#include "portFileio.h"
#include "lineSearch.h"

extern "C"
{
  //NOTE: Note that internally port is using 1 based indexing but
  //externally it is the same as our zero based.. 
  //i.e. c++ a[0] <=> f77 a[1]

  //VARIABLES
  //ndv - number of design variables
  //scale - chosen such that |D[i] * x_i| are all comparable 
  //x - initial guess for optimal vector
  //bounds - upper and lower bounds for x
  //function_ - function to compute cost function
  //gradient_ - function to compute gradient
  //iv - controls internals of algorithm, ie. number of iterations, etc.
  //liv - length of iv, must be > 59+ndv
  //lv - length of v, must be > 71+ndv*(ndv+21)/2
  //v - controls internals of algorithm, ie. length of first step, convergence, etc.
  //uiparm - passed untouched to function_ and gradient_ integer parameters
  //urparm - passed untouched to function_ and gradient_ integer parameters
  extern void dmngb_(int* ndv, double* scale, double* x, double* bounds, void(*function_)(FUNC_ARGS), 
		     void(*gradient_)(GRAD_ARGS), int* iv, int* liv, int* lv, double* v, 
		     int* uiparm, double* urparm);

  //kind - must be set to 2 to alter internal parameters
  extern void divset_(int* kind, int* iv, int* liv, int* lv, double* v);
}

std::string casename;
std::string pathname;
int nprocs;
int locationSwitch;

using namespace std;
int main(int argc, char* argv[])
{
  int i;
  int ndv;

  if(argc != 4){
    std::cerr << "Invalid arguments!!!" << std::endl;
    std::cerr << "USE: " << argv[0] << " <Design filename> <nprocs> <0 - local/ 1 - cluster>" << std::endl;
    return (1);
  }

  casename = argv[1];
  std::string caseMinimal;
  size_t pos = casename.rfind('/');
  if(pos != std::string::npos){
    pathname = casename.substr(0, pos+1);
    caseMinimal = casename.substr(pos);
  }
  else{
    pathname = "./";
  }

  std::stringstream ss;
  ss.clear();
  ss << argv[2];
  ss >> nprocs;
  ss.clear();
  ss << argv[3];
  ss >> locationSwitch;
  
  //check input
  if(locationSwitch != 0 && locationSwitch != 1){
    std::cerr << "WARNING: location switch needs to be 0 - local or 1 - cluster ... received " 
	      << locationSwitch << std::endl;
    return (-1);
  }

  //read first line in design file... not ideal but we need allocations done
  ifstream fin;
  string desfname = casename + ".design";
  fin.open(desfname.c_str());

  if(fin.is_open()){
    fin >> ndv;
    cout << "Number of design variables expected --> " << ndv << endl;
  }
  else{
    std::cerr << "READFILE: Cannot open design file --> " << desfname << std::endl;
    std::cerr << "WARNING: Cannot continue without number of design variables parameter" << std::endl;
    return(1);
  }
  fin.close();

  double* x = new double[ndv];
  double* bounds = new double[2*ndv];
  int* dvType = new int[ndv];
  double* scale = new double[ndv];
  double* v = new double[71 + ndv*(ndv + 19)/2 + 50];
  int lv = 71 + ndv*(ndv + 19)/2 + 50;
  int* iv = new int[59 + ndv];
  int liv = 59 + ndv;
  int* uiparm = new int[ndv+1];         //user int params passed to function_ and gradient_
  double* urparm = new double[ndv+1];   //user real params passed to function_ and gradient_
  
  //initialize iv and v for port usage
  for(i = 0; i < lv; i++){
    v[i] = 0.0;
  }
  for(i = 0; i < liv; i++){
    iv[i] = 0;
  }

  //set iteration limit
  int kind = 2;
  divset_(&kind, iv, &liv, &lv, v);
  //function evaluation limit
  iv[17] = 50;
  //max number of iterations limit
  iv[16] = 50;

  for(i = 0; i < ndv; i++){
    uiparm[i] = 0;
    urparm[i] = 0.0;
    scale[i] = 1.0;
  }
  urparm[ndv] = 0.0;

  std::vector<Param<Real>* > paramList; 
  if(ReadParamFile(paramList, caseMinimal, pathname, false)){
    return (-1);
  }
  for(UInt i = 0; i < paramList.size(); ++i){
    paramList[i]->mode = -1;
  }

  //before we begin we save a copy of the original geometry for all modifications
  SaveOriginalMesh(nprocs);
  SaveOriginalDesignFile();
  
  ReadDesignFile(casename, &ndv, x, bounds, &urparm[0], &urparm[1], dvType);
  //we can intercept the starting values of x here and initialize them
  //but for now we assume the start position is given in the design file
  WriteDesignFile(casename, ndv, x, bounds, urparm[0], &urparm[1], dvType);

  if(0){
    dmngb_(&ndv, scale, x, bounds, function_, gradient_, iv, &liv, &lv, v, uiparm, urparm);
  }
  else{
    lineSearch(ndv, x, bounds, function_, gradient_, uiparm, urparm);
  }

  WriteDesignFile(casename, ndv, x, bounds, urparm[0], &urparm[1], dvType);

  cout << "\nRESULTS:" << endl;
  cout << "=======:" << endl;
  for(i = 0; i < ndv; i++){
    cout << "x[" << i << "] = " << x[i] << endl;
  }

  delete [] x;
  delete [] bounds;
  delete [] scale;
  delete [] v;
  delete [] iv;
  delete [] uiparm;
  delete [] urparm;
  delete [] dvType;

  return (0);
}
