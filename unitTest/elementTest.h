#include <gtest/gtest.h>
#include "elementClass.h"
#include <vector>


class elementTest : public ::testing::Test
{
 public:
  elementTest()
    {
      pointsxyz = new double[4*4];
      //0,0,0
      pointsxyz[0] = 0.0;
      pointsxyz[1] = 0.0;
      pointsxyz[2] = 0.0;
      //1,0,0
      pointsxyz[3] = 1.0;
      pointsxyz[4] = 0.0;
      pointsxyz[5] = 0.0;
      //1,1,0
      pointsxyz[6] = 1.0;
      pointsxyz[7] = 1.0;
      pointsxyz[8] = 0.0;
      //0,1,0
      pointsxyz[9] = 0.0;
      pointsxyz[10] = 1.0;
      pointsxyz[11] = 0.0;
      //0,1,-1
      pointsxyz[12] = 0.0;
      pointsxyz[13] = 1.0;
      pointsxyz[14] = -1.0;
    };
  ~elementTest()
    {
      delete [] pointsxyz;
    };

  //setup points that are spaced on a uniform grid
  double* pointsxyz;
};

TEST_F(elementTest, testTriNormal){
  Element<double>* tri = new Triangle<double>();
  Int listnodes[] = {0,1,2};
  tri->Init(listnodes);
  std::vector<double> normal;
  tri->GetNormal(normal, pointsxyz);
  EXPECT_NEAR(normal[0], 0.0, 1.0e-16);
  EXPECT_NEAR(normal[1], 0.0, 1.0e-16);
  EXPECT_NEAR(normal[2], 1.0, 1.0e-16);
  EXPECT_NEAR(normal[3], 0.5, 1.0e-16);

  EXPECT_NEAR(sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]), 1.0, 1.0e-16);

  delete tri;
}

TEST_F(elementTest, testQuadNormal){
  Element<double>* quad = new Quadrilateral<double>();
  Int listnodes[] = {0,1,2,3};
  quad->Init(listnodes);
  std::vector<double> normal;
  quad->GetNormal(normal, pointsxyz);
  EXPECT_NEAR(normal[0], 0.0, 1.0e-16);
  EXPECT_NEAR(normal[1], 0.0, 1.0e-16);
  EXPECT_NEAR(normal[2], 1.0, 1.0e-16);
  EXPECT_NEAR(normal[3], 1.0, 1.0e-16);

  EXPECT_NEAR(sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]), 1.0, 1.0e-16);
  
  delete quad;

  //now check a deformed quad with a dropped corner
  Element<double>* quadb = new Quadrilateral<double>();
  Int listnodesb[] = {0,1,2,4};
  quadb->Init(listnodesb);
  quad->GetNormal(normal, pointsxyz);
  EXPECT_NEAR(normal[0], -.408, 1.0e-3);
  EXPECT_NEAR(normal[1], 0.408, 1.0e-3);
  EXPECT_NEAR(normal[2], 0.816, 1.0e-3);
  EXPECT_NEAR(normal[3], 1.224, 1.0e-3);

  EXPECT_NEAR(sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]), 1.0, 1.0e-16);
	     
  delete quadb;
}
