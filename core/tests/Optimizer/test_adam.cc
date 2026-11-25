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
#include "Optimist/Optimizer/Adam.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

// Adam is designed for noisy gradients: on deterministic problems it
// converges, but slowly and only to moderate accuracy, hence the looser
// tolerances
using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Trid2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  Function fun;
  Optimist::Optimizer::Adam<Vector> sol;
  sol.task(fun.name());
  sol.learning_rate(0.1);
  sol.tolerance(1.0e-4);
  for (const bool amsgrad : {false, true}) {
    sol.amsgrad_mode(amsgrad);
    SCOPED_TRACE(sol.name());
    expect_solves_all<Problem::OPTIMIZE>(sol, fun, 1.0e-3);
  }
}

TEST(Adam, StepDecay) {
  EllipticParaboloid<Eigen::Vector2d> fun;
  Optimist::Optimizer::Adam<Eigen::Vector2d> sol;
  sol.learning_rate(0.5);
  sol.decay(1.0e-3);
  sol.tolerance(1.0e-4);
  EXPECT_TRUE(solves<Problem::OPTIMIZE>(sol, fun, fun.guess(0), 1.0e-3));
}
