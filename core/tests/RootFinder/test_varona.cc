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
#include "Optimist/RootFinder/Varona.hh"
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

TYPED_TEST(Functions, RootFind) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::Varona<Scalar> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());

  using Method = typename Optimist::RootFinder::Varona<Scalar>::Method;
  for (const auto method : {Method::ORDER_4, Method::ORDER_8, Method::ORDER_16, Method::ORDER_32}) {
    sol.method(method);
  for (const bool damped : {false, true}) {
    sol.damped_mode(damped);
    expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
  }
  }
}

TEST(Varona, ReportsIterationLimit) {
  Sinh<double> fun;
  Optimist::RootFinder::Varona<double> sol;
  sol.max_iterations(2);
  double x_out;
  sol.rootfind(fun, 10.0, x_out);
  EXPECT_FALSE(sol.converged());
  EXPECT_LE(sol.iterations(), 2);
}
