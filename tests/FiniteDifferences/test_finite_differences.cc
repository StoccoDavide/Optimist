/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist/FiniteDifferences.hh"
#include <Eigen/Dense>

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

using Real = double;
constexpr int N = 4;
using Vector = Eigen::Matrix<Real, N, 1>;
using Matrix = Eigen::Matrix<Real, N, N>;

// Structure to hold test functions and their analytical derivatives
struct TestFunction {
  const char* name;
  std::function<bool(const Vector &, Real   &)> scalar_fun;
  std::function<bool(const Vector &, Vector &)> vector_fun;
  std::function<void(const Vector &, Vector &)> grad_analytical;
  std::function<void(const Vector &, Matrix &)> jac_analytical;
  std::function<void(const Vector &, Matrix &)> hes_analytical;
};

// Quadratic scalar and vector functions
TestFunction quadratic = { "Quadratic",

  // Scalar function: f(x) = x(0)^2 + x(1)^2 + x(2)^2 + x(3)^2
  [] (const Vector & x, Real & out) -> bool {
    out = x.squaredNorm();
    return true;
  },

  // Vector function: f(x) = [x_0^2, x_1^2, x_2^2, x_3^2]
  [] (const Vector & x, Vector & out) -> bool {
    for (int i = 0; i < N; ++i) {out(i) = x(i) * x(i);}
    return true;
  },

  // Analytical gradient: 2x
  [] (const Vector & x, Vector & grad) {
    grad = 2.0 * x;
  },

  // Analytical Jacobian: diagonal with 2x_i
  [] (const Vector & x, Matrix & jac) {
    jac.setZero();
    for (int i = 0; i < N; ++i) {jac(i, i) = 2.0 * x(i);};
  },

  // Analytical Hessian: 2*I
  [] (const Vector & /*x*/, Matrix & hes) {
    hes.setZero();
    hes.diagonal().setConstant(2.0);
  }
};

// More complicated scalar and vector functions
TestFunction complicated = {"Complicated",

  // Scalar function: f(x) = x(0)*x(1) + sin(x(2)) + exp(x(3))
  [] (const Vector & x, Real & out) -> bool {
    out = x(0) * x(1) + std::sin(x(2)) + std::exp(x(3));
    return true;
  },

  // Vector function: f(x) = [x(0)*x(1), sin(x(2)), exp(x(3)), x(0)*x(3)]
  [] (const Vector & x, Vector & out) -> bool {
    out(0) = x(0) * x(1);
    out(1) = std::sin(x(2));
    out(2) = std::exp(x(3));
    out(3) = x(0) * x(3);
    return true;
  },

  // Analytical gradient
  [] (const Vector & x, Vector & grad) {
    grad(0) = x(1);
    grad(1) = x(0);
    grad(2) = std::cos(x(2));
    grad(3) = std::exp(x(3));
  },

  // Analytical Jacobian
  [] (const Vector & x, Matrix & jac) {
    jac.setZero();
    jac(0, 0) = x(1);
    jac(0, 1) = x(0);
    jac(1, 2) = std::cos(x(2));
    jac(2, 3) = std::exp(x(3));
    jac(3, 0) = x(3);
    jac(3, 3) = x(0);
  },

  // Analytical Hessian
  [] (const Vector & x, Matrix & hes) {
    hes.setZero();
    hes(0, 1) = 1.0;
    hes(1, 0) = 1.0;
    hes(2, 2) = -std::sin(x(2));
    hes(3, 3) = std::exp(x(3));
  }
};

// Add more test functions here as needed
TestFunction test_functions[] = {quadratic, complicated};

TEST_CASE("FiniteDifferences: Gradient, Jacobian, Hessian", "[finite_differences]") {

  // Tolerance for floating-point comparisons
  Real CBRT_EPSILON{std::pow(std::numeric_limits<Real>::epsilon(), 1.0/3.0)};

  Vector x(Vector::LinSpaced(N, 1.0, N * 1.0));
  for (const auto & tf : test_functions) {
    Vector grad_fd, grad_an;
    Matrix jac_fd, jac_an, hes_fd, hes_an;

    // Gradient
    Optimist::FiniteDifferences::Gradient(x, tf.scalar_fun, grad_fd);
    tf.grad_analytical(x, grad_an);
    REQUIRE((grad_fd - grad_an).norm() < CBRT_EPSILON);

    // Jacobian
    Optimist::FiniteDifferences::Jacobian(x, tf.vector_fun, jac_fd);
    tf.jac_analytical(x, jac_an);
    REQUIRE((jac_fd - jac_an).norm() < CBRT_EPSILON);

    // Hessian
    Optimist::FiniteDifferences::Hessian(x, tf.scalar_fun, hes_fd);
    tf.hes_analytical(x, hes_an);
    REQUIRE((hes_fd - hes_an).norm() < CBRT_EPSILON);
  }
}
