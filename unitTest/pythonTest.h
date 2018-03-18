#include <gtest/gtest.h>
#include <pythonInterface.h>
#include <string>
#include <iostream>
#include <vector>

#ifdef _HAS_PYTHON

class PySetup : public ::testing::Test
{

 protected:
  PySetup(){};
  ~PySetup(){};

  virtual void SetUp()
  {
    pywrap = new PythonWrapper("../unitTest/","testScript","multiply");
  };
  void TearDown()
  {
    delete pywrap;
  };

  PythonWrapper* pywrap;
  
};


class PySetupNumpy : public ::testing::Test
{

 protected:
  PySetupNumpy(){};
  ~PySetupNumpy(){};

  virtual void SetUp()
  {
    pywrap2 = new PythonWrapper("../unitTest/","testScript","arraytest");
  };
  void TearDown()
  {
    delete pywrap2;
  };

  PythonWrapper* pywrap2;
  
};

class PySetupIC : public ::testing::Test
{

 protected:
  PySetupIC(){};
  ~PySetupIC(){};

  virtual void SetUp()
  {
    pywrap2 = new PythonWrapper("../unitTest/","testScript","setInitialConditions");
  };
  void TearDown()
  {
    delete pywrap2;
  };

  PythonWrapper* pywrap2;
  
};

TEST_F(PySetup, testPyMultiply)
{
  int answer = pywrap->CallTwoIntFunction(23, 3);
  EXPECT_EQ(69, answer);
}

TEST_F(PySetupNumpy, testPyArray)
{
  int height = 5;
  std::vector<double> newv(height, 10.0);
  std::vector<double> pReturn = pywrap2->CallDoubleVectorFunction(newv);
  EXPECT_NEAR(0, pReturn[0], 1e-12);
  EXPECT_NEAR(1, pReturn[1], 1e-12);
  EXPECT_NEAR(2, pReturn[2], 1e-12);
  EXPECT_NEAR(3, pReturn[3], 1e-12);
  EXPECT_NEAR(4, pReturn[4], 1e-12);
}

TEST_F(PySetupIC, testInitialConditions)
{
  double* qinf = new double[1];
  double* qret = new double[1];
  double* coords = new double[3];

  qinf[0] = 100.0;
  qret[0] = -9;
  coords[0] = 1.0;
  coords[1] = 2.0;
  coords[2] = 3.0;
  
  pywrap2->SetInitialConditions(qinf, 1, 1, qret, coords);

  EXPECT_NEAR(qret[0], 10.0, 1.0e-12);

  delete [] qinf;
  delete [] qret;
  delete [] coords;
}

#endif
