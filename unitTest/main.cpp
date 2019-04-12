#ifdef _HAS_PYTHON
// we have to include this first b/c python doesn't check for _XOPEN_SOURCE
// and hasn't fixed it yet
#include <Python.h>
#endif

#include <gtest/gtest.h>

//list of all of the headers which define the tests
#include "meshTest.h"
#include "newtonTest.h"
#include "chemTest.h"
#include "compressibleFRTest.h"
#include "gradientTest.h"
#include "pythonTest.h"
#include "pythonInterface.h"
#include "elementTest.h"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  //make HDF quiet
  HDF_TurnOffErrorHandling();
 
  MPI_Init(&argc, &argv);

  Int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Int results = RUN_ALL_TESTS();

  MPI_Finalize();

  return results;
}
