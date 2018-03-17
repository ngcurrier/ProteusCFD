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

#if 1
TEST_F(PySetup, testPyMultiply)
{
  int answer = pywrap->CallTwoIntFunction(23, 3);
  EXPECT_EQ(69, answer);
}
#endif

TEST_F(PySetupNumpy, testPyArray)
{
  int height = 5;
  std::vector<double> newv(height, 10.0);
  std::vector<double> pReturn = pywrap2->CallDoubleVectorFunction(newv);
  EXPECT_NEAR(0, pReturn[0], 10e-12);
  EXPECT_NEAR(1, pReturn[1], 10e-12);
  EXPECT_NEAR(2, pReturn[2], 10e-12);
  EXPECT_NEAR(3, pReturn[3], 10e-12);
  EXPECT_NEAR(4, pReturn[4], 10e-12);
}

#endif
