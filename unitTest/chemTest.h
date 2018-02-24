#include <gtest/gtest.h>
#include <chem.h>
#include <species.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <string>

class ChemSetup : public ::testing::Test
{
 protected:
  ChemSetup(){};
  ~ChemSetup(){};
  std::string homePath;

  virtual void SetUp()
  {
    const char *homedir;
    if ((homedir = getenv("HOME")) == NULL) {
      homedir = getpwuid(getuid())->pw_dir;
    }
    homePath = std::string(homedir);
  };
  void TearDown(){};
};

class ChemTestRead : public ChemSetup
{
 protected:
  ChemTestRead(){};
  ~ChemTestRead(){};

  void SetUp()
  {
    ChemSetup::SetUp();
  };
  void TearDown(){};
};

TEST_F(ChemTestRead, testChemRxnFileRead)
{
  Int isViscous = true;
  std::string chemdb = homePath + "/.proteusCFD/database/chemdb.hdf5";
  ChemModel<double> chem("../unitTest/chemResources/11speciesAir", isViscous, chemdb);
  EXPECT_EQ(11, chem.nspecies);
  EXPECT_EQ(20, chem.nreactions);
  EXPECT_EQ(0, chem.nespecies);

  double rhoi[11];

  //TODO: check speed of sound and mixture properties
  
}

class testSpeciesN2 : public ChemSetup
{
 protected:
  testSpeciesN2(){};
  ~testSpeciesN2(){};
  void SetUp()
  {
    ChemSetup::SetUp();
  };  
  void TearDown(){};
};

TEST_F(testSpeciesN2, testProps)
{
  Species<double> n2;
  std::string chemdb = homePath + "/.proteusCFD/database/chemdb.hdf5";
  n2.Init("N2", true, chemdb);
 
  //1.039 kj/kg.K
  EXPECT_NEAR(1039.6, n2.GetCp(300), 1.0e-1);
  EXPECT_NEAR(1363.6, n2.GetCp(5500), 1.0e-1);
  EXPECT_NEAR(296.8, n2.R, 1.0e-1);
  EXPECT_NEAR(.0280134, n2.MW, 1.0e-7);
  EXPECT_EQ(0.0, n2.charge);
}

class testIdealGasEOS : public ChemSetup
{
 protected:
  testIdealGasEOS(){};
  ~testIdealGasEOS(){};
  void SetUp()
  {
    ChemSetup::SetUp();
  };  
  void TearDown(){};
};

TEST_F(testIdealGasEOS, testCalls)
{
  Species<double> n2;
  std::string chemdb = homePath + "/.proteusCFD/database/chemdb.hdf5";
  n2.Init("N2", true, chemdb);

  IdealGasEOS<double> eos;

  double rho = 1.1384671105; // kg/m^3
  double P = 101325;  // 1 atm
  
  double T = eos.GetT(n2.R, rho, P);
  EXPECT_NEAR(300, T, 0.14);
  double rhocheck = eos.GetRho(n2.R, P, T);
  EXPECT_NEAR(rho, rhocheck, 1.0e-15);
  double pcheck = eos.GetP(n2.R, rho, T);
  EXPECT_NEAR(101325, pcheck, 1.0e-12);
  EXPECT_NEAR(742.9, eos.GetCv(n2.GetCp(300), n2.R, rho, P, T), 1.0e-1);
  EXPECT_NEAR(rho*n2.R, eos.GetRhoR(P, T), 1.0e-12);
}

