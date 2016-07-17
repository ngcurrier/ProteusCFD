#include <gtest/gtest.h>
#include "newton.h"

TEST (NewtonTest, testRootFinder)
{
  //test the one parameter function
  BAR b(32);
  EXPECT_NEAR(3.103569, rootFinder(b, &BAR::f, 2000.0, 1.0e-4), 1.0e-4);

  //test a 3 parameter function
  FOO a;
  Functor<FOO> q(a, &FOO::f, 3.0, 1.0, 32.0);
  EXPECT_NEAR(3.103569, rootFinder(q, &Functor<FOO>::wrapper, 2000.0, 1.0e-4), 1.0e-4);
}
