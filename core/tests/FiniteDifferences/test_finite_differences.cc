/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist/FiniteDifferences.hh"

// Google Test framework
#include <gtest/gtest.h>

using Optimist::TypeTrait;

// Structure to hold test functions and their analytical derivatives
template <typename Vector,
          typename Matrix,
          typename Scalar = typename Vector::Scalar>
  requires TypeTrait<Vector>::IsEigen && TypeTrait<Matrix>::IsEigen &&
           TypeTrait<Scalar>::IsScalar
struct TestFunction {
  using VectorTrait = TypeTrait<Vector>;
  using MatrixTrait = TypeTrait<Matrix>;
  using ScalarTrait = TypeTrait<Scalar>;
  std::function<bool(const Vector &, Scalar &)> scalar_fun;
  std::function<bool(const Vector &, Vector &)> vector_fun;
  std::function<void(const Vector &, Vector &)> grad_analytical;
  std::function<void(const Vector &, Matrix &)> jac_analytical;
  std::function<void(const Vector &, Matrix &)> hes_analytical;
};

// Quadratic scalar and vector functions
template <typename Vector,
          typename Matrix,
          typename Scalar = typename Vector::Scalar>
struct Quadratic : public TestFunction<Vector, Matrix, Scalar> {
  Quadratic() {
    // Scalar function: f(x) = x(0)^2 + x(1)^2 + x(2)^2 + x(3)^2
    this->scalar_fun = [](const Vector &x, Scalar &out) -> bool {
      out = x.squaredNorm();
      return true;
    };

    // Vector function: f(x) = [x_0^2, x_1^2, x_2^2, x_3^2]
    this->vector_fun = [](const Vector &x, Vector &out) -> bool {
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        out.resize(x.size());
      }
      for (Eigen::Index i{0}; i < x.size(); ++i) {
        if constexpr (TypeTrait<Vector>::IsFixed ||
                      TypeTrait<Vector>::IsDynamic) {
          out(i) = x(i) * x(i);
        } else if constexpr (TypeTrait<Vector>::IsSparse) {
          out.coeffRef(i) = x.coeff(i) * x.coeff(i);
        }
      }
      return true;
    };

    // Analytical gradient: 2x
    this->grad_analytical = [](const Vector &x, Vector &grad) {
      grad = 2.0 * x;
    };

    // Analytical Jacobian: diagonal with 2x_i
    this->jac_analytical = [](const Vector &x, Matrix &jac) {
      if constexpr (TypeTrait<Vector>::IsDynamic) {
        jac.resize(x.size(), x.size());
        jac.setZero();
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        jac.resize(x.size(), x.size());
      }
      for (Eigen::Index i{0}; i < x.size(); ++i) {
        if constexpr (TypeTrait<Vector>::IsFixed ||
                      TypeTrait<Vector>::IsDynamic) {
          jac(i, i) = 2.0 * x(i);
        } else if constexpr (TypeTrait<Vector>::IsSparse) {
          jac.coeffRef(i, i) = 2.0 * x.coeff(i);
        }
      };
    };

    // Analytical Hessian: 2*I
    this->hes_analytical = [](const Vector &x, Matrix &hes) {
      if constexpr (TypeTrait<Vector>::IsDynamic) {
        hes.resize(x.size(), x.size());
        hes.setZero();
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        hes.resize(x.size(), x.size());
      }
      for (Eigen::Index i{0}; i < x.size(); ++i) {
        if constexpr (TypeTrait<Matrix>::IsFixed ||
                      TypeTrait<Matrix>::IsDynamic) {
          hes(i, i) = 2.0;
        } else if constexpr (TypeTrait<Matrix>::IsSparse) {
          hes.coeffRef(i, i) = 2.0;
        }
      };
    };
  }
};

// More complex scalar and vector functions
template <typename Vector,
          typename Matrix,
          typename Scalar = typename Vector::Scalar>
struct Complex : public TestFunction<Vector, Matrix, Scalar> {
  Complex() {
    // Scalar function: f(x) = x(0)*x(1) + sin(x(2)) + exp(x(3))
    this->scalar_fun = [](const Vector &x, Scalar &out) -> bool {
      out = x(0) * x(1) + std::sin(x(2)) + std::exp(x(3));
      return true;
    };

    // Vector function: f(x) = [x(0)*x(1), sin(x(2)), exp(x(3)), x(0)*x(3)]
    this->vector_fun = [](const Vector &x, Vector &out) -> bool {
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        out.resize(x.size());
      }
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        out(0) = x(0) * x(1);
        out(1) = std::sin(x(2));
        out(2) = std::exp(x(3));
        out(3) = x(0) * x(3);
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        out.coeffRef(0) = x.coeff(0) * x.coeff(1);
        out.coeffRef(1) = std::sin(x.coeff(2));
        out.coeffRef(2) = std::exp(x.coeff(3));
        out.coeffRef(3) = x.coeff(0) * x.coeff(3);
      }
    };

    // Analytical gradient
    this->grad_analytical = [](const Vector &x, Vector &grad) {
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        grad.resize(x.size());
      }
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        grad(0) = x(1);
        grad(1) = x(0);
        grad(2) = std::cos(x(2));
        grad(3) = std::exp(x(3));
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        grad.coeffRef(0) = x.coeff(1);
        grad.coeffRef(1) = x.coeff(0);
        grad.coeffRef(2) = std::cos(x.coeff(2));
        grad.coeffRef(3) = std::exp(x.coeff(3));
      }
    };

    // Analytical Jacobian
    this->jac_analytical = [](const Vector &x, Matrix &jac) {
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        jac.resize(x.size(), x.size());
      }
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        jac.setZero();
        jac(0, 0) = x(1);
        jac(0, 1) = x(0);
        jac(1, 2) = std::cos(x(2));
        jac(2, 3) = std::exp(x(3));
        jac(3, 0) = x(3);
        jac(3, 3) = x(0);
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        jac.setZero();
        jac.coeffRef(0, 0) = x.coeff(1);
        jac.coeffRef(0, 1) = x.coeff(0);
        jac.coeffRef(1, 2) = std::cos(x.coeff(2));
        jac.coeffRef(2, 3) = std::exp(x.coeff(3));
        jac.coeffRef(3, 0) = x.coeff(3);
        jac.coeffRef(3, 3) = x.coeff(0);
      }
    };

    // Analytical Hessian
    this->hes_analytical = [](const Vector &x, Matrix &hes) {
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        hes.resize(x.size(), x.size());
      }
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        hes.setZero();
        hes(0, 1) = 1.0;
        hes(1, 0) = 1.0;
        hes(2, 2) = -std::sin(x(2));
        hes(3, 3) = std::exp(x(3));
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        hes.setZero();
        hes.coeffRef(0, 1) = 1.0;
        hes.coeffRef(1, 0) = 1.0;
        hes.coeffRef(2, 2) = -std::sin(x.coeff(2));
        hes.coeffRef(3, 3) = std::exp(x.coeff(3));
      }
    };
  }
};

template <typename T>
class FiniteDifferences : public ::testing::Test {};

using TestTypes = ::testing::Types<
    std::pair<Eigen::Vector4d, Eigen::Matrix4d>,
    std::pair<Eigen::VectorXd, Eigen::MatrixXd>,
    std::pair<Eigen::SparseVector<double>, Eigen::SparseMatrix<double>>>;

TYPED_TEST_SUITE(FiniteDifferences, TestTypes);

TYPED_TEST(FiniteDifferences, Quadratic) {
  using Vector = typename TypeParam::first_type;
  using Matrix = typename TypeParam::second_type;
  using Scalar = typename Vector::Scalar;

  // Tolerance for floating-point comparisons
  Scalar CBRT_EPSILON{std::cbrt(std::numeric_limits<Scalar>::epsilon())};

  // Define vector and matrix types
  Vector x;
  if constexpr (TypeTrait<Vector>::IsFixed)
    x = Vector::LinSpaced(Vector::SizeAtCompileTime,
                          1.0,
                          Vector::SizeAtCompileTime);
  else if constexpr (TypeTrait<Vector>::IsDynamic) {
    x.resize(4);
    x = Vector::LinSpaced(4, 1.0, 4.0);
  } else if constexpr (TypeTrait<Vector>::IsSparse) {
    x.resize(4);
    for (size_t i{0}; i < 4; ++i) {
      x.coeffRef(i) = i + 1;
    }
  }

  {
    // Instantiate the quadratic test function
    Quadratic<Vector, Matrix, Scalar> quadratic;

    // Compute gradient
    Vector grad_fd, grad_an;
    bool out_grad{
      Optimist::FiniteDifferences::Gradient(quadratic.scalar_fun, x, grad_fd)};
    EXPECT_TRUE(out_grad);
    quadratic.grad_analytical(x, grad_an);
    EXPECT_TRUE((grad_fd - grad_an).norm() < CBRT_EPSILON);

    // Compute Jacobian
    Matrix jac_fd, jac_an;
    bool out_jac{
      Optimist::FiniteDifferences::Jacobian(quadratic.vector_fun, x, jac_fd)};
    EXPECT_TRUE(out_jac);
    quadratic.jac_analytical(x, jac_an);
    EXPECT_TRUE((jac_fd - jac_an).norm() < CBRT_EPSILON);

    // Compute Hessian
    // Matrix hes_fd, hes_an;
    // bool out_hes{
    //  Optimist::FiniteDifferences::Hessian(quadratic.scalar_fun, x, hes_fd)};
    // EXPECT_TRUE(out_hes);
    // quadratic.hes_analytical(x, hes_an);
    // EXPECT_TRUE((hes_fd - hes_an).norm() < CBRT_EPSILON);
  }

  // Test with the more complex function
  {
    // Instantiate the quadratic test function
    Complex<Vector, Matrix, Scalar> complex;

    // Compute gradient
    Vector grad_fd, grad_an;
    bool out_grad{
      Optimist::FiniteDifferences::Gradient(complex.scalar_fun, x, grad_fd)};
    EXPECT_TRUE(out_grad);
    complex.grad_analytical(x, grad_an);
    EXPECT_TRUE((grad_fd - grad_an).norm() < CBRT_EPSILON);

    // Compute Jacobian
    Matrix jac_fd, jac_an;
    bool out_jac{
      Optimist::FiniteDifferences::Jacobian(complex.vector_fun, x, jac_fd)};
    EXPECT_TRUE(out_jac);
    complex.jac_analytical(x, jac_an);
    EXPECT_TRUE((jac_fd - jac_an).norm() < CBRT_EPSILON);

    // Compute Hessian
    // Matrix hes_fd, hes_an;
    // bool out_hes{
    //   Optimist::FiniteDifferences::Hessian(complex.scalar_fun, x, hes_fd)};
    // EXPECT_TRUE(out_hes);
    // complex.hes_analytical(x, hes_an);
    // EXPECT_TRUE((hes_fd - hes_an).norm() < CBRT_EPSILON);
  }
}
