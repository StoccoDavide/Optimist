/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Checks on the test functions themselves: analytic derivatives against
// central finite differences, and known solutions being actual roots or
// stationary points.

// Optimist library
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;
using Optimist::Utils::from_dense;
using Optimist::Utils::to_dense;

namespace {

  // Relative agreement between an analytic and a finite-difference quantity
  template <typename A, typename B>
  testing::AssertionResult agree(const A &analytic,
                                 const B &numeric,
                                 const double tol,
                                 const std::string &what) {
    const double err{static_cast<double>((analytic - numeric).norm())};
    const double scale{1.0 + static_cast<double>(analytic.norm())};
    if (err <= tol * scale) {
      return testing::AssertionSuccess();
    }
    return testing::AssertionFailure()
           << what << " mismatch: |analytic - fd| = " << err
           << " > " << tol * scale << "\nanalytic:\n"
           << analytic << "\nfinite differences:\n"
           << numeric;
  }

  constexpr double FD_TOLERANCE{1.0e-5};

  // Fourth-order central difference of a vector-valued map along the j-th
  // coordinate (Richardson extrapolation of the classic central difference)
  template <typename Map>
  Eigen::VectorXd richardson(Map &&map,
                             const Eigen::VectorXd &x,
                             const Integer j,
                             const double h) {
    auto shifted = [&](const double t) {
      Eigen::VectorXd y(x);
      y(j) += t;
      return Eigen::VectorXd(map(y));
    };
    return (8 * (shifted(h) - shifted(-h)) - (shifted(2 * h) - shifted(-2 * h))) /
           (12 * h);
  }

}  // namespace

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Scalar functions f : R -> R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

using ScalarTestTypes = ConcatT<ScalarVariants<Linear>,
                                ScalarVariants<Quadratic>,
                                ScalarVariants<Cos>,
                                ScalarVariants<Cosh>,
                                ScalarVariants<Sin>,
                                ScalarVariants<Sinh>>;

OPTIMIST_TYPED_TEST_SUITE(ScalarFunctions, ScalarTestTypes);

TYPED_TEST(ScalarFunctions, DerivativesMatchFiniteDifferences) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;
  Function fun;
  for (const Scalar x0 : fun.guesses()) {
    const double x{static_cast<double>(x0)};
    const double h{std::cbrt(static_cast<double>(Function::EPSILON)) *
                   std::max(1.0, std::abs(x))};
    Scalar f_p, f_m, d_p, d_m, d, dd;
    ASSERT_TRUE(fun.evaluate(x0 + h, f_p));
    ASSERT_TRUE(fun.evaluate(x0 - h, f_m));
    ASSERT_TRUE(fun.first_derivative(x0, d));
    ASSERT_TRUE(fun.first_derivative(x0 + h, d_p));
    ASSERT_TRUE(fun.first_derivative(x0 - h, d_m));
    ASSERT_TRUE(fun.second_derivative(x0, dd));
    const double tol{std::cbrt(static_cast<double>(Function::EPSILON))};
    EXPECT_NEAR(d, (f_p - f_m) / (2 * h), tol * (1 + std::abs(d)));
    EXPECT_NEAR(dd, (d_p - d_m) / (2 * h), tol * (1 + std::abs(dd)));
  }
}

TYPED_TEST(ScalarFunctions, SolutionsAreRootsOrStationaryPoints) {
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;
  Function fun;
  ASSERT_FALSE(fun.solutions().empty());
  ASSERT_FALSE(fun.guesses().empty());
  const Scalar tol{100 * Function::SQRT_EPSILON};
  for (const Scalar s : fun.solutions()) {
    Scalar f, d;
    ASSERT_TRUE(fun.evaluate(s, f));
    ASSERT_TRUE(fun.first_derivative(s, d));
    EXPECT_TRUE(std::abs(f) < tol || std::abs(d) < tol)
        << fun.name() << " at " << s << ": f = " << f << ", f' = " << d;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Square nonlinear systems F : R^n -> R^n
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

using SystemTestTypes = ConcatT<DoubleVectorVariants<Linear1, 1>,
                                DoubleVectorVariants<Booth, 2>,
                                DoubleVectorVariants<Rosenbrock2, 2>,
                                DoubleVectorVariants<Rosenbrock4, 4>,
                                DoubleVectorVariants<Rosenbrock8, 8>,
                                DoubleVectorVariants<PowellSingular, 4>>;

OPTIMIST_TYPED_TEST_SUITE(SystemFunctions, SystemTestTypes);

TYPED_TEST(SystemFunctions, JacobianMatchesFiniteDifferences) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Dense    = Eigen::VectorXd;
  Function fun;
  for (const Vector &x0 : fun.guesses()) {
    const Dense x(to_dense(x0));
    const Integer n{static_cast<Integer>(x.size())};
    typename Function::FirstDerivative J;
    ASSERT_TRUE(fun.first_derivative(x0, J));
    Eigen::MatrixXd J_fd(n, n);
    for (Integer j{0}; j < n; ++j) {
      const double h{1.0e-5};
      auto F = [&fun](const Dense &y) {
        Vector out;
        EXPECT_TRUE(fun.evaluate(from_dense<Vector>(y), out));
        return Dense(to_dense(out));
      };
      J_fd.col(j) = richardson(F, x, j, h);
    }
    EXPECT_TRUE(agree(Eigen::MatrixXd(J), J_fd, FD_TOLERANCE, "Jacobian"));
  }
}

TYPED_TEST(SystemFunctions, SolutionsAreRoots) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  Function fun;
  ASSERT_FALSE(fun.solutions().empty());
  ASSERT_FALSE(fun.guesses().empty());
  for (const Vector &s : fun.solutions()) {
    Vector F;
    ASSERT_TRUE(fun.evaluate(s, F));
    EXPECT_LT(F.norm(), 1.0e-10) << fun.name() << " at " << to_string(s);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Cost functions f : R^n -> R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

using CostTestTypes = ConcatT<DoubleVectorVariants<EllipticParaboloid, 2>,
                              DoubleVectorVariants<Schaffer2, 2>,
                              DoubleVectorVariants<Himmelblau, 2>,
                              DoubleVectorVariants<Beale, 2>,
                              DoubleVectorVariants<Trid2, 2>,
                              DoubleVectorVariants<Trid6, 6>,
                              DoubleVectorVariants<RosenbrockCost2, 2>,
                              DoubleVectorVariants<RosenbrockCost4, 4>,
                              DoubleVectorVariants<Rastrigin2, 2>,
                              DoubleVectorVariants<Ackley2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(CostFunctions, CostTestTypes);

TYPED_TEST(CostFunctions, DerivativesMatchFiniteDifferences) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Dense    = Eigen::VectorXd;
  Function fun;
  for (const Vector &x0 : fun.guesses()) {
    const Dense x(to_dense(x0));
    const Integer n{static_cast<Integer>(x.size())};
    typename Function::FirstDerivative g;
    typename Function::SecondDerivative H;
    ASSERT_TRUE(fun.first_derivative(x0, g));
    ASSERT_TRUE(fun.second_derivative(x0, H));
    Dense g_fd(n);
    Eigen::MatrixXd H_fd(n, n);
    for (Integer j{0}; j < n; ++j) {
      const double h{1.0e-5};
      auto f = [&fun](const Dense &y) {
        double out;
        EXPECT_TRUE(fun.evaluate(from_dense<Vector>(y), out));
        return Dense::Constant(1, out).eval();
      };
      auto grad = [&fun](const Dense &y) {
        typename Function::FirstDerivative out;
        EXPECT_TRUE(fun.first_derivative(from_dense<Vector>(y), out));
        return Dense(to_dense(out));
      };
      g_fd(j)     = richardson(f, x, j, h)(0);
      H_fd.col(j) = richardson(grad, x, j, h);
    }
    EXPECT_TRUE(agree(to_dense(g), g_fd, FD_TOLERANCE, "gradient"))
        << fun.name() << " at " << to_string(x0);
    EXPECT_TRUE(agree(Eigen::MatrixXd(H), H_fd, FD_TOLERANCE, "Hessian"))
        << fun.name() << " at " << to_string(x0);
  }
}

TYPED_TEST(CostFunctions, SolutionsAreStationaryPoints) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  Function fun;
  ASSERT_FALSE(fun.solutions().empty());
  ASSERT_FALSE(fun.guesses().empty());
  for (const Vector &s : fun.solutions()) {
    typename Function::FirstDerivative g;
    ASSERT_TRUE(fun.first_derivative(s, g));
    EXPECT_LT(g.norm(), 1.0e-9) << fun.name() << " at " << to_string(s);
  }
}

TEST(CostFunctions, TridMinimumValue) {
  Trid6<Eigen::VectorXd> fun;
  double f;
  ASSERT_TRUE(fun.evaluate(fun.solution(0), f));
  EXPECT_NEAR(f, Trid6<Eigen::VectorXd>::minimum(), 1.0e-12);
  EXPECT_LT(f, 0.0);
}

TEST(SystemFunctions, BrownIsOverdetermined) {
  using In  = Eigen::Vector2d;
  using Out = Eigen::Vector3d;
  Brown<In, Out> fun;
  Out F;
  ASSERT_TRUE(fun.evaluate(fun.solution(0), F));
  EXPECT_LT(F.norm(), 1.0e-10);

  // Jacobian (3 x 2) and component Hessians (three 2 x 2 matrices)
  const Eigen::VectorXd x(to_dense(fun.guess(0)));
  Brown<In, Out>::FirstDerivative J;
  Brown<In, Out>::SecondDerivative H;
  ASSERT_TRUE(fun.first_derivative(fun.guess(0), J));
  ASSERT_TRUE(fun.second_derivative(fun.guess(0), H));
  ASSERT_EQ(H.size(), 3u);
  Eigen::MatrixXd J_fd(3, 2);
  for (Integer j{0}; j < 2; ++j) {
    auto map = [&fun](const Eigen::VectorXd &y) {
      Out out;
      EXPECT_TRUE(fun.evaluate(from_dense<In>(y), out));
      return Eigen::VectorXd(out);
    };
    J_fd.col(j) = richardson(map, x, j, 1.0e-5);
  }
  EXPECT_TRUE(agree(Eigen::MatrixXd(J), J_fd, FD_TOLERANCE, "Jacobian"));
  EXPECT_EQ(H[2](0, 1), 1.0);
  EXPECT_EQ(H[2](1, 0), 1.0);
  EXPECT_EQ(H[0].norm(), 0.0);
}
