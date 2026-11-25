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
#include "Optimist/Optimizer/GoldenSection.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<ScalarVariants<Quadratic>,
                          ScalarVariants<Cos>,
                          ScalarVariants<Sin>,
                          ScalarVariants<Cosh>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::Optimizer::GoldenSection<Scalar> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  for (const Scalar x_ini : fun.guesses()) {
    Scalar x_out;
    EXPECT_TRUE(minimizes(sol, fun, x_ini, &x_out));
    Scalar dd;
    ASSERT_TRUE(fun.second_derivative(x_out, dd));
    EXPECT_GT(dd, 0) << "not a minimum";
  }
}

TEST(GoldenSection, BoundedSearch) {
  // On [0, 6] cos(x) is unimodal with minimum at pi
  Cos<double> fun;
  Optimist::Optimizer::GoldenSection<double> sol;
  sol.bounds(0.0, 6.0);
  sol.tolerance(solver_tolerance<double>());
  double x_out;
  sol.optimize(fun, 1.0, x_out);
  EXPECT_TRUE(sol.converged());
  EXPECT_NEAR(x_out, M_PI, value_tolerance<double>());
}

TEST(GoldenSection, UnboundedBelowIsNotConverged) {
  // A line has no minimum: the bracketing phase must give up
  Linear<double> fun;
  Optimist::Optimizer::GoldenSection<double> sol;
  sol.max_function_evaluations(10000);
  double x_out;
  EXPECT_FALSE(sol.optimize(fun, 0.0, x_out));
}
