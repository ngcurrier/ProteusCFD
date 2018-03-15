#include <gtest/gtest.h>
#include <pythonInterface.h>

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

TEST_F(PySetup, testPyMultiply)
{
  int answer = pywrap->CallTwoIntFunction(23, 3);
  EXPECT_EQ(69, answer);
}

#endif
