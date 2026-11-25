/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// One-dimensional root finders handed an optimization problem look for the
// stationary points of phi = f^2/2, i.e., the roots of f*f'

// Optimist library
#include "Optimist/RootFinder/NewtonRaphson.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<ScalarVariants<Linear>,
                          ScalarVariants<Quadratic>,
                          ScalarVariants<Cos>,
                          ScalarVariants<Sin>,
                          ScalarVariants<Sinh>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, StationaryPointsOfLeastSquares) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::NewtonRaphson<Scalar> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  // Only the damped iteration is tested: the undamped one cycles between 0
  // and -1 on Quadratic (a property of Newton's method on f*f')
  sol.enable_damped_mode();
  expect_solves_all<Problem::OPTIMIZE>(sol, fun, default_tolerance<Scalar>());
}
