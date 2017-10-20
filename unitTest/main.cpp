#include <gtest/gtest.h>

//list of all of the headers which define the tests
#include "meshTest.h"
#include "newtonTest.h"
#include "chemTest.h"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  Int rank, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Int results = RUN_ALL_TESTS();

  MPI_Finalize();

  return results;
}
