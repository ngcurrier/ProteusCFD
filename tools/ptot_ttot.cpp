#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 6){
    cout << "USAGE: " << argv[0] << " <Pstatic> <Tstatic> <gamma> <Uref> <Uinf>" << endl; 
    return 1;
  }
  
  stringstream ss (stringstream::in | stringstream::out);

  double p;
  double t;
  double gamma;
  double uref;
  double uinf;

  cout.precision(12);
  cout.setf(ios_base::scientific);

  ss.clear();
  ss << argv[1];
  ss >> p;

  ss.clear();
  ss << argv[2];
  ss >> t;

  ss.clear();
  ss << argv[3];
  ss >> gamma;

  ss.clear();
  ss << argv[4];
  ss >> uref;

  ss.clear();
  ss << argv[5];
  ss >> uinf;

  cout << "\nGiven: " << endl;
  cout << "=======:" << endl;
  cout << "Pstatic: " << p << endl;
  cout << "Tstatic: " << t << endl;
  cout << "gamma: " << gamma << endl;
  cout << "Uref: " << uref << endl;
  cout << "Uinf: " << uinf << endl;


  double Ma = uinf/uref;

  cout << "\nComputed: " << endl;
  cout << "=========" << endl;
  cout << "Mach: " << Ma << endl;

  double gm1 = gamma -1.0;
  double power = gamma/gm1;
  double factor = (1.0+(gm1/2.0)*Ma*Ma);
  double pt = p*pow(factor,power);
  double tt = t*factor;

  cout << "Ptotal: " << pt << endl;
  cout << "Ttotal: " << tt << endl;

  return 0;
}
