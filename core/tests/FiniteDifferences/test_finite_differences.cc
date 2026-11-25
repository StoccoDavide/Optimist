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

#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

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
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        jac.resize(x.size(), x.size());
      }
      jac.setZero();
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
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        hes.resize(x.size(), x.size());
      }
      hes.setZero();
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
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        out = x(0) * x(1) + std::sin(x(2)) + std::exp(x(3));
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        out = x.coeff(0) * x.coeff(1) + std::sin(x.coeff(2)) +
              std::exp(x.coeff(3));
      }
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
      return true;
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
      jac.setZero();
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        jac(0, 0) = x(1);
        jac(0, 1) = x(0);
        jac(1, 2) = std::cos(x(2));
        jac(2, 3) = std::exp(x(3));
        jac(3, 0) = x(3);
        jac(3, 3) = x(0);
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
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
      hes.setZero();
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsFixed) {
        hes(0, 1) = 1.0;
        hes(1, 0) = 1.0;
        hes(2, 2) = -std::sin(x(2));
        hes(3, 3) = std::exp(x(3));
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
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
    EXPECT_LT((grad_fd - grad_an).norm(), CBRT_EPSILON);

    // Compute Jacobian
    Matrix jac_fd, jac_an;
    bool out_jac{
      Optimist::FiniteDifferences::Jacobian(quadratic.vector_fun, x, jac_fd)};
    EXPECT_TRUE(out_jac);
    quadratic.jac_analytical(x, jac_an);
    EXPECT_LT((jac_fd - jac_an).norm(), CBRT_EPSILON);

    // Compute Hessian
    Matrix hes_fd, hes_an;
    bool out_hes{
      Optimist::FiniteDifferences::Hessian(quadratic.scalar_fun, x, hes_fd)};
    EXPECT_TRUE(out_hes);
    quadratic.hes_analytical(x, hes_an);
    EXPECT_LT((hes_fd - hes_an).norm(), CBRT_EPSILON);
  }

  {
    // Instantiate the complex test function
    Complex<Vector, Matrix, Scalar> complex;

    // Compute gradient
    Vector grad_fd, grad_an;
    bool out_grad{
      Optimist::FiniteDifferences::Gradient(complex.scalar_fun, x, grad_fd)};
    EXPECT_TRUE(out_grad);
    complex.grad_analytical(x, grad_an);
    EXPECT_LT((grad_fd - grad_an).norm(), CBRT_EPSILON);

    // Compute Jacobian
    Matrix jac_fd, jac_an;
    bool out_jac{
      Optimist::FiniteDifferences::Jacobian(complex.vector_fun, x, jac_fd)};
    EXPECT_TRUE(out_jac);
    complex.jac_analytical(x, jac_an);
    EXPECT_LT((jac_fd - jac_an).norm(), CBRT_EPSILON);

    // Compute Hessian
    Matrix hes_fd, hes_an;
    bool out_hes{
      Optimist::FiniteDifferences::Hessian(complex.scalar_fun, x, hes_fd)};
    EXPECT_TRUE(out_hes);
    complex.hes_analytical(x, hes_an);
    EXPECT_LT((hes_fd - hes_an).norm(), CBRT_EPSILON);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Finite differences against the analytic derivatives of the test set
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

using namespace OptimistTest;
using Optimist::Utils::to_dense;
namespace TS = Optimist::TestSet;

namespace {

  // Relative error between two dense objects
  template <typename A, typename B>
  double relative_error(const A &a, const B &b) {
    using Dense = Eigen::Matrix<typename A::Scalar, Eigen::Dynamic,
                                Eigen::Dynamic>;
    const Eigen::MatrixXd da{Dense(a).template cast<double>()},
        db{Dense(b).template cast<double>()};
    return (da - db).norm() / (1.0 + da.norm());
  }

  // Tolerance of the centered differences: O(eps^(2/3)), plus a safety
  // factor for the curvature of the test functions
  template <typename Scalar>
  double fd_tolerance() {
    return 1.0e3 * std::pow(static_cast<double>(
                                 std::numeric_limits<Scalar>::epsilon()),
                             2.0 / 3.0);
  }

}  // namespace

using CostTypes = ConcatT<VectorVariants<TS::EllipticParaboloid, 2>,
                          VectorVariants<TS::Himmelblau, 2>,
                          VectorVariants<TS::Beale, 2>,
                          VectorVariants<TS::Trid6, 6>,
                          VectorVariants<TS::RosenbrockCost4, 4>,
                          DoubleVectorVariants<TS::Rastrigin2, 2>>;

OPTIMIST_TYPED_TEST_SUITE(CostFunctions, CostTypes);

TYPED_TEST(CostFunctions, GradientAndHessian) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;
  using Matrix   = typename Function::SecondDerivative;

  Function fun;
  auto f = [&fun](const Vector &x, Scalar &out) {
    return fun.evaluate(x, out);
  };
  for (const Vector &x : fun.guesses()) {
    typename Function::FirstDerivative g;
    Matrix H;
    ASSERT_TRUE(fun.first_derivative(x, g));
    ASSERT_TRUE(fun.second_derivative(x, H));
    Vector g_fd;
    Matrix H_fd;
    ASSERT_TRUE(Optimist::FiniteDifferences::Gradient(f, x, g_fd));
    ASSERT_TRUE(Optimist::FiniteDifferences::Hessian(f, x, H_fd));
    EXPECT_LT(relative_error(to_dense(g), to_dense(g_fd)),
              fd_tolerance<Scalar>())
        << fun.name() << " at " << to_string(x);
    // Second differences of function values: O(sqrt(eps)) accuracy
    EXPECT_LT(relative_error(H, H_fd),
              100 * std::sqrt(static_cast<double>(Function::EPSILON)))
        << fun.name() << " at " << to_string(x);
  }
}

using SystemTypes = ConcatT<VectorVariants<TS::Booth, 2>,
                            VectorVariants<TS::Rosenbrock4, 4>,
                            VectorVariants<TS::PowellSingular, 4>>;

OPTIMIST_TYPED_TEST_SUITE(SystemFunctions, SystemTypes);

TYPED_TEST(SystemFunctions, Jacobian) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;
  using Matrix   = typename Function::FirstDerivative;

  Function fun;
  auto F = [&fun](const Vector &x, Vector &out) {
    return fun.evaluate(x, out);
  };
  for (const Vector &x : fun.guesses()) {
    Matrix J, J_fd;
    ASSERT_TRUE(fun.first_derivative(x, J));
    ASSERT_TRUE(Optimist::FiniteDifferences::Jacobian(F, x, J_fd));
    EXPECT_LT(relative_error(J, J_fd), fd_tolerance<Scalar>())
        << fun.name() << " at " << to_string(x);
  }
}

TEST(FiniteDifferences, NonFiniteFunctionIsReported) {
  auto f = [](const Eigen::Vector2d &x, double &out) {
    out = std::log(x(0)) + x(1);  // -inf at x_0 = 0
    return true;
  };
  Eigen::Vector2d g;
  EXPECT_FALSE(Optimist::FiniteDifferences::Gradient(
      f, Eigen::Vector2d(0.0, 1.0), g));
}

TEST(FiniteDifferences, OneSidedNearSingularity) {
  // log(x) is undefined for x < 0: close to the boundary the left point is
  // not finite and the derivative must come from the right side only
  auto f = [](const Eigen::Vector<double, 1> &x, double &out) {
    out = std::log(x(0));
    return std::isfinite(out);
  };
  const Eigen::Vector<double, 1> x(1.0e-9);
  Eigen::Vector<double, 1> g;
  ASSERT_TRUE(Optimist::FiniteDifferences::Gradient(f, x, g));
  EXPECT_NEAR(g(0) * x(0), 1.0, 1.0e-2);
}
