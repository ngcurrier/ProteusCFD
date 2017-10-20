#include <gtest/gtest.h>
#include <chem.h>

class ChemTestRead : public testing::Test
{
 protected:
  ChemTestRead(){};
  ~ChemTestRead(){};

  void SetUp(){};
  void TearDown(){};

};

TEST_F(ChemTestRead, testChemFileRead)
{
  //read the file and check goodness
}
