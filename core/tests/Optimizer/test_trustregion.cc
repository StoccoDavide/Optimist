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
#include "Optimist/Optimizer/TrustRegion.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Himmelblau, 2>,
                          VectorVariants<Beale, 2>,
                          VectorVariants<Trid2, 2>,
                          VectorVariants<Trid6, 6>,
                          VectorVariants<RosenbrockCost2, 2>,
                          DoubleVectorVariants<RosenbrockCost4, 4>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;
  using Solver   = Optimist::Optimizer::TrustRegion<Vector>;
  using Method   = typename Solver::Method;

  Function fun;
  Solver sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  sol.max_iterations(500);
  for (const auto method : {Method::STEIHAUG, Method::DOGLEG}) {
    sol.method(method);
    SCOPED_TRACE(sol.name());
    expect_solves_all<Problem::OPTIMIZE>(sol, fun, default_tolerance<Scalar>());
  }
}

// Vector-valued functions are minimized in the least-squares sense, using
// the full Hessian J'J + sum_i F_i H_i
using SystemTypes = ConcatT<VectorVariants<Booth, 2>,
                            VectorVariants<Rosenbrock2, 2>,
                            DoubleVectorVariants<Rosenbrock4, 4>>;

OPTIMIST_TYPED_TEST_SUITE(Systems, SystemTypes);

TYPED_TEST(Systems, LeastSquares) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::Optimizer::TrustRegion<Vector> sol;
  sol.tolerance(std::pow(Function::EPSILON, Scalar(0.75)));
  sol.max_iterations(500);
  expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
}

TEST(TrustRegion, NegativeCurvatureAtStart) {
  // At (0, 0) the Hessian of Himmelblau's function is indefinite: Newton's
  // step is not a descent direction, the trust region keeps the iteration
  // safe and the dogleg falls back to Steihaug-Toint
  Himmelblau<Eigen::Vector2d> fun;
  for (const auto method :
       {Optimist::Optimizer::TrustRegion<Eigen::Vector2d>::Method::STEIHAUG,
        Optimist::Optimizer::TrustRegion<Eigen::Vector2d>::Method::DOGLEG}) {
    Optimist::Optimizer::TrustRegion<Eigen::Vector2d> sol;
    sol.method(method);
    EXPECT_TRUE(solves<Problem::OPTIMIZE>(sol, fun, Eigen::Vector2d(0.0, 0.0),
                                          default_tolerance<double>()))
        << sol.name();
  }
}
