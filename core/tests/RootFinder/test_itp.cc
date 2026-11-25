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
#include "Optimist/RootFinder/ITP.hh"
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

TYPED_TEST(Functions, RootFindInBracket) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::ITP<Scalar> sol;
  sol.bounds(0.0, static_cast<Scalar>(M_PI));
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
}

TYPED_TEST(Functions, RootFindInShiftedBracket) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  // Asymmetric bracket around a known root
  Function fun;
  Optimist::RootFinder::ITP<Scalar> sol;
  sol.tolerance(solver_tolerance<Scalar>());
  for (const Scalar root : fun.solutions()) {
    Scalar f;
    ASSERT_TRUE(fun.evaluate(root, f));
    if (std::abs(f) > Function::SQRT_EPSILON) {
      continue;  // Stationary point, not a root
    }
    sol.bounds(root - Scalar(0.3), root + Scalar(0.7));
    Scalar x_out;
    sol.rootfind(fun, root, x_out);
    EXPECT_TRUE(sol.converged()) << fun.name() << " around " << root;
    EXPECT_NEAR(x_out, root, default_tolerance<Scalar>())
        << fun.name() << " around " << root;
  }
}

TEST(ITP, NoSignChangeIsNotConverged) {
  // cosh(x) > 0 everywhere: nothing is bracketed
  Cosh<double> fun;
  Optimist::RootFinder::ITP<double> sol;
  sol.bounds(-1.0, 2.0);
  double x_out;
  EXPECT_FALSE(sol.rootfind(fun, 0.0, x_out));
  EXPECT_FALSE(sol.converged());
}
