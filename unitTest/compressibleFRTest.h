#include <solutionSpace.h>
#include <param.h>
#include <string>
#include <vector>

TEST(compressibleFR, newtonTest)
{
  SolutionSpace<double>* space = NULL;
  std::string pathname = "../unitTest/chemResources/";
  std::string casestring = "flatplate";

  std::vector<Param<Real>* > paramList; 
  ReadParamFile(paramList, casestring, pathname);

  CompressibleFREqnSet<double>* eqnset = new CompressibleFREqnSet<double>(space, paramList[0]);

  double rPress = paramList[0]->ref_pressure;
  double rTemp = paramList[0]->ref_temperature;
  double rRho = paramList[0]->ref_density;
  
  double Pgoal = 101200/rPress;
  
  // 79% N2 and 21% O2
  double rho = 1.161/rRho; //kg/m^3
  double rhoi[5] = {0.21*rho, 0.0, 0.79*rho, 0.0, 0.0};
  
  double Treturn = eqnset->NewtonFindTGivenP(rhoi, Pgoal, 150.0/rTemp);
  Treturn *= rTemp;
  EXPECT_NEAR(301.57, Treturn, 0.01);

  Treturn = eqnset->NewtonFindTGivenP(rhoi, Pgoal/2.0, 150.0/rTemp);
  Treturn *= rTemp;
  EXPECT_NEAR(301.57/2.0, Treturn, 0.01);

  delete eqnset;
}
