/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Behaviour shared by all solvers through SolverBase: bounds, counters,
// limits, verbose output, problem semantics and the lambda interface.

// Optimist library
#include "Optimist/Optimizer/BFGS.hh"
#include "Optimist/Optimizer/NelderMead.hh"
#include "Optimist/RootFinder/Halley.hh"
#include "Optimist/RootFinder/Newton.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;
using Optimist::Optimizer::BFGS;
using Optimist::Optimizer::NelderMead;
using Optimist::RootFinder::Halley;
using Optimist::RootFinder::Newton;
using Optimist::RootFinder::NewtonRaphson;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Bounds
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TEST(Bounds, DefaultsAreInfinite) {
  NewtonRaphson<double> scalar;
  EXPECT_EQ(scalar.lower_bound(), -std::numeric_limits<double>::infinity());
  EXPECT_EQ(scalar.upper_bound(), +std::numeric_limits<double>::infinity());
  EXPECT_FALSE(scalar.has_finite_bounds());

  Newton<Eigen::Vector2d> fixed;
  EXPECT_TRUE((fixed.lower_bound().array() < 0).all());
  EXPECT_FALSE(fixed.lower_bound().allFinite());
  EXPECT_FALSE(fixed.has_finite_bounds());

  // Dynamic and sparse bounds are empty until set
  Newton<Eigen::VectorXd> dynamic;
  EXPECT_EQ(dynamic.lower_bound().size(), 0);
  EXPECT_FALSE(dynamic.has_finite_bounds());
  Newton<Eigen::SparseVector<double>> sparse;
  EXPECT_FALSE(sparse.has_finite_bounds());
}

TEST(Bounds, ValidBoundsAreAccepted) {
  NewtonRaphson<double> scalar;
  EXPECT_NO_THROW(scalar.bounds(-1.0, 2.0));
  EXPECT_TRUE(scalar.has_finite_bounds());

  Newton<Eigen::Vector2d> fixed;
  EXPECT_NO_THROW(fixed.bounds(Eigen::Vector2d(-1, -2), Eigen::Vector2d(1, 2)));
  EXPECT_TRUE(fixed.has_finite_bounds());
  EXPECT_NO_THROW(fixed.upper_bound(Eigen::Vector2d(3, 3)));
  EXPECT_NO_THROW(fixed.lower_bound(Eigen::Vector2d(0, 0)));

  Newton<Eigen::SparseVector<double>> sparse;
  EXPECT_NO_THROW(sparse.bounds(make<Eigen::SparseVector<double>>({-1, -1}),
                                make<Eigen::SparseVector<double>>({1, 1})));
  EXPECT_TRUE(sparse.has_finite_bounds());
}

TEST(Bounds, InvalidBoundsAreRejected) {
  NewtonRaphson<double> scalar;
  EXPECT_THROW(scalar.bounds(2.0, -1.0), std::runtime_error);
  EXPECT_THROW(scalar.bounds(1.0, 1.0), std::runtime_error);

  Newton<Eigen::Vector2d> fixed;
  // Only one component is inverted
  EXPECT_THROW(fixed.bounds(Eigen::Vector2d(-1, 3), Eigen::Vector2d(1, 2)),
               std::runtime_error);
  fixed.bounds(Eigen::Vector2d(-1, -1), Eigen::Vector2d(1, 1));
  EXPECT_THROW(fixed.lower_bound(Eigen::Vector2d(0, 2)), std::runtime_error);
  EXPECT_THROW(fixed.upper_bound(Eigen::Vector2d(0, -2)), std::runtime_error);

  Newton<Eigen::VectorXd> dynamic;
  EXPECT_THROW(dynamic.bounds(Eigen::VectorXd::Zero(2), Eigen::VectorXd::Ones(3)),
               std::runtime_error);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Counters and limits
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TEST(Counters, ResetBetweenSolves) {
  Rosenbrock2<Eigen::Vector2d> fun;
  Newton<Eigen::Vector2d> sol;
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  const Integer iterations{sol.iterations()};
  const Integer evaluations{sol.function_evaluations()};
  const Integer jacobians{sol.jacobian_evaluations()};
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  EXPECT_EQ(sol.iterations(), iterations);
  EXPECT_EQ(sol.function_evaluations(), evaluations);
  EXPECT_EQ(sol.jacobian_evaluations(), jacobians);
  EXPECT_GT(jacobians, 0);
}

TEST(Counters, HessianCountsSecondDerivatives) {
  Cos<double> fun;
  Halley<double> sol;
  double x_out;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  EXPECT_GT(sol.hessian_evaluations(), 0);
  EXPECT_EQ(sol.hessian_evaluations(), sol.jacobian_evaluations());
}

TEST(Limits, IterationLimitIsReported) {
  RosenbrockCost2<Eigen::Vector2d> fun;
  BFGS<Eigen::Vector2d> sol;
  sol.max_iterations(3);
  EXPECT_EQ(sol.max_iterations(), 3);
  Eigen::Vector2d x_out;
  EXPECT_FALSE(sol.optimize(fun, fun.guess(0), x_out));
  EXPECT_LE(sol.iterations(), 3);
}

TEST(Limits, FunctionEvaluationLimitThrowsInDerivativeBasedSolvers) {
  Rosenbrock2<Eigen::Vector2d> fun;
  Newton<Eigen::Vector2d> sol;
  sol.max_function_evaluations(1);
  Eigen::Vector2d x_out;
  EXPECT_THROW(sol.rootfind(fun, fun.guess(0), x_out), std::runtime_error);
}

TEST(Limits, FunctionEvaluationLimitStopsDerivativeFreeSolvers) {
  RosenbrockCost2<Eigen::Vector2d> fun;
  NelderMead<Eigen::Vector2d> sol;
  sol.max_function_evaluations(15);
  Eigen::Vector2d x_out;
  // The simplex method relies on the solver bookkeeping: the limit is
  // enforced by an exception raised at the 16-th evaluation
  EXPECT_THROW(sol.optimize(fun, fun.guess(0), x_out), std::runtime_error);
  EXPECT_LE(sol.function_evaluations(), 15);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Verbose output
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TEST(Verbose, TableIsWrittenToTheChosenStream) {
  Rosenbrock2<Eigen::Vector2d> fun;
  Newton<Eigen::Vector2d> sol;
  std::ostringstream os;
  sol.ostream(os);
  sol.task("MyTask");
  sol.enable_verbose_mode();
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  const std::string out{os.str()};
  EXPECT_NE(out.find("Newton"), std::string::npos);
  EXPECT_NE(out.find("MyTask"), std::string::npos);
  EXPECT_NE(out.find("CONVERGED"), std::string::npos);
  EXPECT_EQ(out.find("NOT CONVERGED"), std::string::npos);
}

TEST(Verbose, SilentByDefault) {
  Rosenbrock2<Eigen::Vector2d> fun;
  Newton<Eigen::Vector2d> sol;
  std::ostringstream os;
  sol.ostream(os);
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  EXPECT_TRUE(os.str().empty());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Problem semantics
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TEST(Semantics, OptimizersMinimizeCostFunctionsDirectly) {
  // min Trid2 = -2 at (2, 2). Minimizing f^2/2 would stop anywhere on the
  // ellipse f = 0
  Trid2<Eigen::Vector2d> fun;
  BFGS<Eigen::Vector2d> sol;
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.optimize(fun, Eigen::Vector2d(0, 0), x_out));
  EXPECT_NEAR(x_out(0), 2.0, 1.0e-6);
  EXPECT_NEAR(x_out(1), 2.0, 1.0e-6);
}

TEST(Semantics, OptimizersRootFindThroughLeastSquares) {
  Booth<Eigen::Vector2d> fun;
  BFGS<Eigen::Vector2d> sol;
  sol.tolerance(1.0e-12);
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), x_out));
  EXPECT_TRUE(fun.is_solution(x_out, 1.0e-6));
}

TEST(Semantics, RootFindersOptimizeThroughStationaryPoints) {
  // Newton's method on grad f = 0 with Jacobian = Hessian of f: one step on a
  // quadratic
  Trid2<Eigen::Vector2d> fun;
  Newton<Eigen::Vector2d> sol;
  sol.disable_damped_mode();
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.optimize(fun, Eigen::Vector2d(0, 0), x_out));
  EXPECT_NEAR(x_out(0), 2.0, 1.0e-10);
  EXPECT_NEAR(x_out(1), 2.0, 1.0e-10);
  EXPECT_LE(sol.iterations(), 3);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Lambda interface
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TEST(Lambdas, SolveWithUserLambdas) {
  // f(x) = (x_0 - 1)^2 + 10 (x_1 + 2)^2
  auto f = [](const Eigen::Vector2d &x, double &out) {
    out = (x(0) - 1) * (x(0) - 1) + 10 * (x(1) + 2) * (x(1) + 2);
    return true;
  };
  auto g = [](const Eigen::Vector2d &x, Eigen::Vector2d &out) {
    out << 2 * (x(0) - 1), 20 * (x(1) + 2);
    return true;
  };
  BFGS<Eigen::Vector2d> sol;
  Eigen::Vector2d x_out;
  ASSERT_TRUE(sol.solve(f, g, Eigen::Vector2d(5, 5), x_out));
  EXPECT_NEAR(x_out(0), 1.0, 1.0e-6);
  EXPECT_NEAR(x_out(1), -2.0, 1.0e-6);

  NelderMead<Eigen::Vector2d> simplex;
  simplex.max_iterations(1000);
  ASSERT_TRUE(simplex.solve(f, Eigen::Vector2d(5, 5), x_out));
  EXPECT_NEAR(x_out(0), 1.0, 1.0e-6);
  EXPECT_NEAR(x_out(1), -2.0, 1.0e-6);
}
