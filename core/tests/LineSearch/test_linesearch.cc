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
#include "Optimist/LineSearch.hh"

// Test utilities
#include "TestUtils.hh"

namespace LS = Optimist::LineSearch;

namespace {

  // phi(a) = (a - 2)^2 - 1, phi(0) = 3, phi'(0) = -4, minimum at a = 2
  struct Parabola {
    int evaluations{0};
    bool operator()(const double a, double &v, double *d) {
      ++evaluations;
      v = (a - 2) * (a - 2) - 1;
      if (d != nullptr) {
        *d = 2 * (a - 2);
      }
      return true;
    }
  };

  // phi(a) = Rosenbrock's function along the steepest descent direction at
  // (-1.2, 1): strongly curved, the unit step is far too long
  struct RosenbrockLine {
    Eigen::Vector2d x{-1.2, 1.0}, d;
    RosenbrockLine() {
      d = -gradient(x);
    }
    static double value(const Eigen::Vector2d &z) {
      return 100 * std::pow(z(1) - z(0) * z(0), 2) + std::pow(1 - z(0), 2);
    }
    static Eigen::Vector2d gradient(const Eigen::Vector2d &z) {
      return Eigen::Vector2d(-400 * z(0) * (z(1) - z(0) * z(0)) - 2 * (1 - z(0)),
                             200 * (z(1) - z(0) * z(0)));
    }
    bool operator()(const double a, double &v, double *dv) const {
      const Eigen::Vector2d z{x + a * d};
      v = value(z);
      if (dv != nullptr) {
        *dv = gradient(z).dot(d);
      }
      return true;
    }
    double phi_0() const {
      return value(x);
    }
    double dphi_0() const {
      return gradient(x).dot(d);
    }
  };

}  // namespace

TEST(LineSearch, ArmijoSufficientDecrease) {
  Parabola phi;
  LS::Settings<double> s;
  double alpha{10.0}, value;
  ASSERT_TRUE(LS::armijo(phi, 3.0, -4.0, alpha, value, s));
  EXPECT_LE(value, 3.0 + s.c1 * alpha * -4.0);
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 10.0);
}

TEST(LineSearch, ArmijoAcceptsTheFirstGoodStep) {
  Parabola phi;
  double alpha{1.0}, value;
  ASSERT_TRUE(LS::armijo(phi, 3.0, -4.0, alpha, value));
  EXPECT_DOUBLE_EQ(alpha, 1.0);
  EXPECT_EQ(phi.evaluations, 1);
}

TEST(LineSearch, GoldsteinConditions) {
  for (const double alpha_0 : {0.01, 1.0, 100.0}) {
    Parabola phi;
    LS::Settings<double> s;
    double alpha{alpha_0}, value;
    ASSERT_TRUE(LS::goldstein(phi, 3.0, -4.0, alpha, value, s)) << alpha_0;
    EXPECT_LE(value, 3.0 + s.goldstein * alpha * -4.0);
    EXPECT_GE(value, 3.0 + (1 - s.goldstein) * alpha * -4.0);
  }
}

TEST(LineSearch, StrongWolfeConditions) {
  for (const double alpha_0 : {0.01, 1.0, 100.0}) {
    Parabola phi;
    LS::Settings<double> s;
    s.c2 = 0.1;
    double alpha{alpha_0}, value, slope;
    ASSERT_TRUE(LS::strong_wolfe(phi, 3.0, -4.0, alpha, value, slope, s))
        << alpha_0;
    EXPECT_LE(value, 3.0 + s.c1 * alpha * -4.0);
    EXPECT_LE(std::abs(slope), s.c2 * 4.0);
  }
}

TEST(LineSearch, MoreThuenteConditions) {
  for (const double alpha_0 : {0.01, 1.0, 100.0}) {
    Parabola phi;
    LS::Settings<double> s;
    s.c2 = 0.1;
    double alpha{alpha_0}, value, slope;
    ASSERT_TRUE(LS::more_thuente(phi, 3.0, -4.0, alpha, value, slope, s))
        << alpha_0;
    EXPECT_LE(value, 3.0 + s.c1 * alpha * -4.0);
    EXPECT_LE(std::abs(slope), s.c2 * 4.0);
  }
}

TEST(LineSearch, MoreThuenteOnCurvedValley) {
  RosenbrockLine phi;
  LS::Settings<double> s;
  double alpha{1.0}, value, slope;
  ASSERT_TRUE(LS::more_thuente(phi, phi.phi_0(), phi.dphi_0(), alpha, value,
                               slope, s));
  EXPECT_LE(value, phi.phi_0() + s.c1 * alpha * phi.dphi_0());
  EXPECT_LE(std::abs(slope), -s.c2 * phi.dphi_0());
}

TEST(LineSearch, MoreThuenteIsCheaperThanStrongWolfeOnAQuadratic) {
  // From a far trial step both searches must zoom: the safeguarded cubic and
  // quadratic steps of More-Thuente find the exact minimizer of a quadratic
  Parabola p_1, p_2;
  LS::Settings<double> s;
  s.c2 = 0.1;
  double a_1{50.0}, a_2{50.0}, v, d;
  ASSERT_TRUE(LS::more_thuente(p_1, 3.0, -4.0, a_1, v, d, s));
  ASSERT_TRUE(LS::strong_wolfe(p_2, 3.0, -4.0, a_2, v, d, s));
  EXPECT_LE(p_1.evaluations, p_2.evaluations + 1);
}

TEST(LineSearch, StrongWolfeOnCurvedValley) {
  RosenbrockLine phi;
  LS::Settings<double> s;
  double alpha{1.0}, value, slope;
  ASSERT_TRUE(LS::strong_wolfe(phi, phi.phi_0(), phi.dphi_0(), alpha, value,
                               slope, s));
  EXPECT_LE(value, phi.phi_0() + s.c1 * alpha * phi.dphi_0());
  EXPECT_LE(std::abs(slope), -s.c2 * phi.dphi_0());
}

TEST(LineSearch, RejectAscentDirections) {
  Parabola phi;
  double alpha{1.0}, value, slope;
  EXPECT_FALSE(LS::armijo(phi, 3.0, +4.0, alpha, value));
  EXPECT_FALSE(LS::goldstein(phi, 3.0, +4.0, alpha, value));
  EXPECT_FALSE(LS::strong_wolfe(phi, 3.0, +4.0, alpha, value, slope));
  EXPECT_FALSE(LS::more_thuente(phi, 3.0, +4.0, alpha, value, slope));
  EXPECT_EQ(phi.evaluations, 0);
}

TEST(LineSearch, CubicInterpolationIsExactOnCubics) {
  // phi(a) = a^3 - 3a, local minimum at a = 1
  auto p  = [](const double a) { return a * a * a - 3 * a; };
  auto dp = [](const double a) { return 3 * a * a - 3; };
  EXPECT_NEAR(LS::cubic_minimizer(0.0, p(0.0), dp(0.0), 3.0, p(3.0), dp(3.0)),
              1.0, 1.0e-12);
  EXPECT_NEAR(LS::cubic_minimizer(3.0, p(3.0), dp(3.0), 0.0, p(0.0), dp(0.0)),
              1.0, 1.0e-12);
}

TEST(LineSearch, CubicInterpolationIsSafeguarded) {
  // The minimizer of the interpolant is outside [0, 1]: the result stays in
  // the inner part of the interval
  const double a{LS::cubic_minimizer(0.0, 0.0, -1.0, 1.0, -2.0, -1.0)};
  EXPECT_GE(a, 0.1);
  EXPECT_LE(a, 0.9);
}
