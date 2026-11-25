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
#include "Optimist/Optimizer/BFGS.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

// Trid<6> is tested in double precision only: with |f| ~ 50 a line search
// based on function values cannot locate the minimizer better than about
// sqrt(eps*|f|/lambda_min) ~ 1e-2 in single precision
using TestTypes = ConcatT<VectorVariants<EllipticParaboloid, 2>,
                          VectorVariants<Himmelblau, 2>,
                          VectorVariants<Beale, 2>,
                          VectorVariants<Trid2, 2>,
                          DoubleVectorVariants<Trid6, 6>,
                          VectorVariants<RosenbrockCost2, 2>,
                          DoubleVectorVariants<RosenbrockCost4, 4>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, LocalMinimum) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;
  using Solver   = Optimist::Optimizer::BFGS<Vector>;
  using Method   = typename Solver::Method;
  using Search   = Optimist::LineSearch::Method;

  Function fun;
  Solver sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  sol.max_iterations(500);
  sol.max_function_evaluations(10000);
  sol.max_gradient_evaluations(10000);

  for (const auto method : {Method::BFGS, Method::LBFGS}) {
    for (const auto search :
         {Search::STRONG_WOLFE, Search::MORE_THUENTE, Search::ARMIJO,
          Search::GOLDSTEIN}) {
      SCOPED_TRACE(Optimist::LineSearch::name(search));
      sol.method(method);
      sol.line_search(search);
      expect_solves_all<Problem::OPTIMIZE>(sol, fun,
                                           default_tolerance<Scalar>());
    }
  }
}

// Vector-valued functions are minimized in the least-squares sense
using SystemTypes = ConcatT<VectorVariants<Booth, 2>,
                            VectorVariants<Rosenbrock2, 2>,
                            DoubleVectorVariants<Rosenbrock4, 4>>;

OPTIMIST_TYPED_TEST_SUITE(Systems, SystemTypes);

TYPED_TEST(Systems, LeastSquares) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::Optimizer::BFGS<Vector> sol;
  // The gradient of |F|^2/2 is J'F, which vanishes linearly at the root: the
  // tolerance on it bounds the error on x up to the conditioning of J
  sol.tolerance(std::pow(Function::EPSILON, Scalar(0.75)));
  sol.max_iterations(500);
  expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
}

TEST(BFGS, MinimizesTheCostDirectly) {
  // The minimum of Trid is negative: minimizing f^2/2 instead of f would end
  // on the zero level set, far from the minimizer
  Trid6<Eigen::VectorXd> fun;
  Optimist::Optimizer::BFGS<Eigen::VectorXd> sol;
  Eigen::VectorXd x_out;
  ASSERT_TRUE(sol.optimize(fun, fun.guess(0), x_out));
  double f;
  ASSERT_TRUE(fun.evaluate(x_out, f));
  EXPECT_NEAR(f, Trid6<Eigen::VectorXd>::minimum(), 1.0e-8);
}

TEST(BFGS, LimitedMemoryOnLargeProblem) {
  // 50-dimensional Trid with only 5 stored pairs. The cost is evaluated with
  // heavy cancellation (terms ~1e7 summing to ~-2e4), so its rounding noise is
  // far above the default line search slack: the relative slack is raised to
  // the value used by Hager and Zhang
  Trid<Eigen::VectorXd, 50> fun;
  Optimist::Optimizer::BFGS<Eigen::VectorXd> sol;
  sol.method(Optimist::Optimizer::BFGS<Eigen::VectorXd>::Method::LBFGS);
  sol.memory(5);
  sol.line_search_settings().f_slack = 1.0e-6;
  sol.tolerance(1.0e-10);
  sol.max_iterations(2000);
  sol.max_function_evaluations(10000);
  sol.max_gradient_evaluations(10000);
  // The solution entries are up to 650 and the Hessian condition number is
  // ~1e3: check the relative error
  Eigen::VectorXd x_out;
  ASSERT_TRUE(sol.optimize(fun, fun.guess(0), x_out));
  EXPECT_LT((x_out - fun.solution(0)).norm() / fun.solution(0).norm(), 1.0e-10);
}
