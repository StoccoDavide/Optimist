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
#include "Optimist/Optimizer/PatternSearch.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

// Coordinate search stalls in the curved valley of Rosenbrock's function in
// single precision, so that case is tested in double precision only
using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Himmelblau, 2>,
                          VectorVariants<Trid2, 2>,
                          DoubleVectorVariants<Trid6, 6>,
                          DoubleVectorVariants<RosenbrockCost2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::Optimizer::PatternSearch<Vector> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  sol.max_iterations(5000);
  sol.max_function_evaluations(50000);
  expect_solves_all<Problem::OPTIMIZE>(sol, fun, default_tolerance<Scalar>());
}

TEST(PatternSearch, BudgetExhaustedIsNotConverged) {
  RosenbrockCost2<Eigen::Vector2d> fun;
  Optimist::Optimizer::PatternSearch<Eigen::Vector2d> sol;
  sol.max_function_evaluations(20);
  Eigen::Vector2d x_out;
  EXPECT_NO_THROW(sol.optimize(fun, fun.guess(0), x_out));
  EXPECT_FALSE(sol.converged());
  EXPECT_LE(sol.function_evaluations(), 20);
}
