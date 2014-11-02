#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
  int i;
  
  ifstream exact;
  ifstream approx;

  string fileexact;
  string fileapprox;

  stringstream ss;

  int size;
  string line;

  double* esol;
  double* asol;
  double* atime;
  double* etime;

  if(argc != 3){
     cerr << "USAGE: " << argv[0] << " filename_exact filename_approx" << endl;
     return (-1);
  }
  ss.clear();
  ss.str() = "";
  ss << argv[1];
  ss >> fileexact;
  ss.clear();
  ss.str() = "";
  ss << argv[2];
  ss >> fileapprox;
  
  //open and read solution file
  exact.open(fileexact.c_str());

  //count number of lines of data
  while(!exact.eof()){
    getline(exact, line);
    size++;
  }
  //subtract blank line and header
  size -= 2;
  cout << "Npts in file: " << size << endl;

  esol = new double[size];
  asol = new double[size];
  etime = new double[size];
  atime = new double[size];

  //seek beginning of file and clear iostream
  exact.clear();
  exact.seekg(0, ios::beg);
  //clear header line
  getline(exact, line);

  //open and read approximate solution file
  approx.open(fileapprox.c_str());
  //clear header line
  getline(approx, line);

  for(i = 0; i < size; i++){
    approx >> atime[i];
    approx >> asol[i];
    //dump two values;
    approx >> line;
    approx >> line;
    exact >> etime[i];
    exact >> esol[i];
    //dump two values;
    exact >> line;
    exact >> line;
  }

  //check that the read was in sync
  for(i = 0; i < size; i++){
    if(atime[i] != etime[i]){
      cerr << "TIME OUT OF SYNC " << atime[i] << " " << etime[i] << endl;
    }
  }

  //compute l2 norm of the error
  double l2 = 0.0;
  double err;
  for(i = 0; i < size; i++){
    err = abs(esol[i] - asol[i]);
    l2 += err*err;
  }
  //  l2 = sqrt(l2)/(double)size;
  l2 = sqrt(l2);

  cout << "L2 ERROR: " << l2 << endl;


  approx.close();
  exact.close();

  delete [] esol;
  delete [] asol;
  delete [] etime;
  delete [] atime;

  return 0;
}

