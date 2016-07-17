#include <gtest/gtest.h>

//list of all of the headers which define the tests
#include "meshTest.h"
#include "newtonTest.h"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
