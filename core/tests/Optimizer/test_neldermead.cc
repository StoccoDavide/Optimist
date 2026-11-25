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
#include "Optimist/Optimizer/NelderMead.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

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
  using Solver   = Optimist::Optimizer::NelderMead<Vector>;
  using Method   = typename Solver::Method;

  Function fun;
  Solver sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  sol.max_iterations(5000);
  sol.max_function_evaluations(20000);
  for (const auto method : {Method::STANDARD, Method::SPENDLEY}) {
    for (const bool adaptive : {false, true}) {
      SCOPED_TRACE(testing::Message() << "method " << static_cast<int>(method)
                                      << ", adaptive " << adaptive);
      sol.method(method);
      sol.adaptive_mode(adaptive);
      expect_minimizes_all(sol, fun);
    }
  }
}

TEST(NelderMead, OneDimensional) {
  // A simplex in one dimension is a segment
  using Vector = Eigen::Vector<double, 1>;
  Rastrigin<Vector, 1> fun;
  Optimist::Optimizer::NelderMead<Vector> sol;
  sol.delta(0.1);
  Vector x_out;
  sol.optimize(fun, Vector(0.2), x_out);
  EXPECT_TRUE(sol.converged());
  EXPECT_NEAR(x_out(0), 0.0, 1.0e-6);
}
