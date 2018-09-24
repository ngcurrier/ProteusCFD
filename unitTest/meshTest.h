#include <gtest/gtest.h>
#include "mesh.h"
#include "parallel.h"

class MeshTestGmsh : public testing::Test
{
 protected:
 MeshTestGmsh():
  ierr(0)
    { };
  ~MeshTestGmsh(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

class MeshTestGmsh4 : public testing::Test
{
 protected:
 MeshTestGmsh4():
  ierr(0)
    { };
  ~MeshTestGmsh4(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

#ifdef _HAS_CGNS
class MeshTestCGNS : public testing::Test
{
 protected:
 MeshTestCGNS():
  ierr(0)
    { };
  ~MeshTestCGNS(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};
#endif

class MeshTestUgrid : public testing::Test
{
 protected:
 MeshTestUgrid():
  ierr(0)
    { };
  ~MeshTestUgrid(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

class MeshTestCrunch : public testing::Test
{
 protected:
 MeshTestCrunch():
  ierr(0)
    { };
  ~MeshTestCrunch(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

class MeshTestSU2 : public testing::Test
{
 protected:
 MeshTestSU2():
  ierr(0)
    { };
  ~MeshTestSU2(){};

  void SetUp(){};
  void TearDown(){};

  int ierr;
  Mesh<Real> m;
};

TEST_F(MeshTestGmsh, testGmshRead)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadGMSH2_Ascii("../unitTest/meshResources/cube.msh"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(14, m.GetNumNodes());
  EXPECT_EQ(48, m.GetNumElem());

  m.BuildMaps();
  m.CalcMetrics();
  EXPECT_NEAR(1.0, m.GetVolumeTotal(), 1.0e-7);
}

TEST_F(MeshTestGmsh4, testGmshRead)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadGMSH4_Ascii("../unitTest/meshResources/Domain.msh"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(56601,m.GetNumNodes());
  EXPECT_EQ(327552, m.GetNumElem());

  m.BuildMaps();
  m.CalcMetrics();
  EXPECT_NEAR(0.00201639, m.GetVolumeTotal(), 1.0e-7);
}

#ifdef _HAS_CGNS
TEST_F(MeshTestCGNS, testCGNSRead)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadCGNS("../unitTest/meshResources/CGNS/demo.cgns"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(21632,m.GetNumNodes());
  EXPECT_EQ(19296, m.GetNumElem());

  //m.BuildMaps();
  //m.CalcMetrics();
  //EXPECT_NEAR(0.00201639, m.GetVolumeTotal(), 1.0e-7);
}
#endif

TEST_F(MeshTestUgrid, testUgridRead)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadUGRID_Ascii("../unitTest/meshResources/nhex.ugrid"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(4276, m.GetNumNodes());
  EXPECT_EQ(10244, m.GetNumElem());
  
  m.BuildMaps();
  m.CalcMetrics();
  EXPECT_NEAR(0.00990099, m.GetVolumeTotal(), 1.0e-7);
}

TEST_F(MeshTestCrunch, testCrunchRead)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadCRUNCH_Ascii("../unitTest/meshResources/naca0012_two_plane.crunch"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(5188, m.GetNumNodes());
  EXPECT_EQ(14724, m.GetNumElem());

  m.BuildMaps();
  m.CalcMetrics();
  EXPECT_NEAR(14.92148669, m.GetVolumeTotal(), 1.0e-7);
}

TEST_F(MeshTestSU2, testSU2Read)
{
  PObj<Real> pobj;
  EXPECT_EQ(0, m.ReadSU2_Ascii("../unitTest/meshResources/mesh_ONERAM6_inv.su2"));
  m.SetParallelPointer(&pobj);

  EXPECT_EQ(108396, m.GetNumNodes());
  EXPECT_EQ(621508, m.GetNumElem());

  m.BuildMaps();
  m.CalcMetrics();
  EXPECT_NEAR(4442.0775193109512, m.GetVolumeTotal(), 1.0e-7);
}
