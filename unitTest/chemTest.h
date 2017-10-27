#include <gtest/gtest.h>
#include <chem.h>
#include <species.h>

class ChemTestRead : public testing::Test
{
 protected:
  ChemTestRead(){};
  ~ChemTestRead(){};

  void SetUp(){};
  void TearDown(){};

};

TEST_F(ChemTestRead, testChemRxnFileRead)
{
  Int isViscous = true;
  std::string chemdb = "/usr/local/database/chemdb.hdf5";
  ChemModel<double> chem("../unitTest/chemResources/11speciesAir", isViscous, chemdb);
  EXPECT_EQ(11, chem.nspecies);
  EXPECT_EQ(20, chem.nreactions);
  EXPECT_EQ(1, chem.nespecies);
}

TEST(testSpeciesN2, cp)
{
  Species<double> n2;
  n2.Init("N2", true, "/usr/local/database/chemdb.hdf5");

  //1.039 kj/kg.K
  EXPECT_NEAR(1039.6, n2.GetCp(300), 1.0e-1);
}
