/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist/RootFinder/Greenstadt.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<VectorVariants<Linear1, 1>,
                          VectorVariants<Booth, 2>,
                          VectorVariants<Rosenbrock2, 2>,
                          VectorVariants<Rosenbrock4, 4>,
                          VectorVariants<Rosenbrock6, 6>,
                          VectorVariants<Rosenbrock8, 8>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, RootFind) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::Greenstadt<Vector> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  using Method = typename Optimist::RootFinder::Greenstadt<Vector>::Method;
  for (const auto method : {Method::ONE, Method::TWO}) {
    sol.method(method);
  for (const bool damped : {false, true}) {
    sol.damped_mode(damped);
    expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
  }
  }
}
