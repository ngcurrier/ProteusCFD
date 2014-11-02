#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include "portFileio.h"

int ReadDesignFile(std::string casename, int* ndv, double* x, double* bounds, double* f, 
		   double* grad, int* dvType)
{
  int i, j, id;
  int ndv_orig;
  int trashid;
  double temp;
  int* map;
  std::ifstream fin;
  std::string filename = casename + ".design";

  fin.open(filename.c_str());

  if(fin.is_open()){
    fin >> (*ndv);
    ndv_orig = (*ndv);
    map = new int[*ndv];
    j = 0;
    for(i = 0; i < ndv_orig; i++){
      fin >> id;
      if(id != i){
	std::cerr << "WARNING: line id " << i << " missing from file" << std::endl;
	return(99);
      }
      fin >> dvType[j];
      fin >> x[j];
      fin >> bounds[j*2 + 0];
      fin >> bounds[j*2 + 1];
      if(bounds[j*2 + 1] < bounds[j*2 + 0]){
	//eliminate design variable from list if bounds are backwards
	//user control tool...
	(*ndv)--;
	map[i] = 0;
      }
      else{
	//increment counter to load next value
	j++;
	map[i] = 1;
      }
    }
    fin >> (*f);
    j = 0;
    for(i = 0; i < ndv_orig; i++){
      fin >> trashid;
      fin >> temp;
      if(map[i]){
	grad[j] = temp;
	j++;
      }
    }
  }
  else{
    std::cerr << "DESIGN READFILE: Cannot open design file --> " << filename << std::endl;
    return(1);
  }

  delete [] map; 

  fin.close();

  return(0);
} 

int WriteDesignFile(std::string casename, int ndv, double* x, double* bounds, double f, 
		    double* grad, int* dvType)
{
  int i;
  std::ofstream fout;
  std::string filename = casename + ".design";
  
  fout.open(filename.c_str());

  if(fout.is_open()){
    fout << ndv << std::endl;
    

    fout.setf(std::ios::scientific);
    fout.precision(15);

    for(i = 0; i < ndv; i++){
      //write out current value of design variable and lower and upper bounds on variable
      //using fortran style array indexing here for port compatibility
      fout << i << " " << dvType[i] << " " << x[i] << " " 
	   << bounds[0 + 2*i] << " " << bounds[1 + 2*i] << std::endl; 
    }

    //write out current value of cost function
    fout << f << std::endl;

    //write out derivatives of cost function w.r.t. design variable
    for(i = 0; i < ndv; i++){
      fout << i << " " << grad[i] << std::endl;
    }

    fout << "#Format is:\n";
    fout << "#number of design variables\n";
    fout << "#beta 0 <designType> <value> <lower bound> <upper bound>\n";
    fout << "#beta 1 <designType> <value> <lower bound> <upper bound>\n";
    fout << "#...\n";
    fout << "#Objective function value\n";
    fout << "#Grad beta 0\n";
    fout << "#Grad beta 1\n";
    fout << "#...\n";

  }
  else{
    std::cerr << "DESIGN WRITEFILE: Cannot open design file --> " << filename << std::endl;
    return(1);
  }

  fout.close();
  return(0);
}

int WriteGradToDesignFile(std::string casename, double* grad)
{
  int err = 0;
  int ndv = GetNdvDesignFile(casename);
  double* x = new double[ndv];
  double* bounds = new double[ndv*2];
  double f;
  double* trashgrad = new double[ndv];
  int* dvType = new int[ndv];

  err = ReadDesignFile(casename, &ndv, x, bounds, &f, trashgrad, dvType);
  err = WriteDesignFile(casename, ndv, x, bounds, f, grad, dvType);

  delete [] x;
  delete [] bounds;
  delete [] trashgrad;
  delete [] dvType;

  return err;
}

int WriteFunctionToDesignFile(std::string casename, double f)
{
  int err = 0;
  int ndv = GetNdvDesignFile(casename);
  double* x = new double[ndv];
  double* bounds = new double[ndv*2];
  double f_old;
  double* grad = new double[ndv];
  int* dvType = new int[ndv];

  err = ReadDesignFile(casename, &ndv, x, bounds, &f_old, grad, dvType);

  err = WriteDesignFile(casename, ndv, x, bounds, f, grad, dvType);

  delete [] x;
  delete [] bounds;
  delete [] grad;
  delete [] dvType;

  return err;
}

int GetNdvDesignFile(std::string casename)
{
  int ndv = -1;
  std::ifstream fin;
  std::string filename = casename + ".design";

  fin.open(filename.c_str());
  if(fin.is_open()){
    fin >> ndv;
  }
  else{
    std::cerr << "DESIGN WRITEFILE: Cannot open design file --> " << filename << std::endl;
    return(-1);
  }
  fin.close();
  return ndv;
}


void ClearLog(std::string casename)
{
  std::string logname = casename + ".designlog";
  std::ofstream logout;
  logout.open(logname.c_str());
  logout.close();
}

void WriteLog(std::string casename, int ndv, double* x, double* grad, double f)
{
  std::string logname = casename + ".designlog";
  std::ofstream logout;
  logout.open(logname.c_str(), std::fstream::app);

  //write current state
  for(int i = 0; i < ndv; i++){
    logout << std::setprecision(16) << x[i] << " ";
  }

  //write gradients
  for(int i = 0; i < ndv; i++){
    logout << std::setprecision(16) << grad[i] << " ";
  }

  //write function evaluation
  logout << std::setprecision(16) << f << std::endl;

  logout.close();
}
