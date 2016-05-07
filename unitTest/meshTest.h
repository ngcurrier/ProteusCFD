#include <gtest/gtest.h>
#include "mesh.h"

class MeshTest : public testing::Test
{
 protected:
 MeshTest():
  ierr(0)
    { };
  ~MeshTest(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

TEST_F(MeshTest, testGmshRead)
{
  EXPECT_EQ(0, m.ReadGMSH_Ascii("test.msh"));
}
