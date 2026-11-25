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
#include "Optimist/RootFinder/LevenbergMarquardt.hh"
#include "Optimist/TestSet.hh"

// Test utilities
#include "TestUtils.hh"

using namespace Optimist::TestSet;
using namespace OptimistTest;

using TestTypes = ConcatT<VectorVariants<Linear1, 1>,
                          VectorVariants<Booth, 2>,
                          VectorVariants<Rosenbrock2, 2>,
                          VectorVariants<Rosenbrock4, 4>,
                          VectorVariants<Rosenbrock6, 6>,
                          VectorVariants<Rosenbrock8, 8>>;

OPTIMIST_TYPED_TEST_SUITE(Functions, TestTypes);

TYPED_TEST(Functions, RootFind) {
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;
  using Scalar   = typename Function::Scalar;

  Function fun;
  Optimist::RootFinder::LevenbergMarquardt<Vector> sol;
  sol.task(fun.name());
  sol.tolerance(solver_tolerance<Scalar>());
  for (const bool scaling : {false, true}) {
    sol.scaling_mode(scaling);
  for (const bool damped : {true}) {
    sol.damped_mode(damped);
    expect_solves_all<Problem::ROOTFIND>(sol, fun, default_tolerance<Scalar>());
  }
  }
}

TEST(LevenbergMarquardt, SingularRoot) {
  // The Jacobian is singular at the root: the convergence is only linear,
  // so the solution is accurate to roughly the square root of the residual
  // tolerance
  PowellSingular<Eigen::Vector4d> fun;
  Optimist::RootFinder::LevenbergMarquardt<Eigen::Vector4d> sol;
  sol.max_iterations(500);
  sol.max_function_evaluations(2000);
  sol.tolerance(1.0e-10);
  Eigen::Vector4d x_out;
  EXPECT_TRUE(solves<Problem::ROOTFIND>(sol, fun, fun.guess(0), 1.0e-3, &x_out));
}

TEST(LevenbergMarquardt, SingularJacobianAtStart) {
  // Rosenbrock's Jacobian is singular on x_1 = 0 only through a zero
  // pivot; Newton's method breaks down while the damped normal equations
  // remain solvable
  Rosenbrock2<Eigen::Vector2d> fun;
  Optimist::RootFinder::LevenbergMarquardt<Eigen::Vector2d> sol;
  sol.max_iterations(200);
  sol.tolerance(solver_tolerance<double>());
  EXPECT_TRUE(solves<Problem::ROOTFIND>(sol, fun, Eigen::Vector2d(0.0, 0.0), default_tolerance<double>()));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Overdetermined systems
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Input, typename Output>
struct BrownCase {
  using In  = Input;
  using Out = Output;
};

using BrownTypes = testing::Types<
    BrownCase<Eigen::Vector<double, 2>, Eigen::Vector<double, 3>>,
    BrownCase<Eigen::VectorXd, Eigen::VectorXd>,
    BrownCase<Eigen::SparseVector<double>, Eigen::SparseVector<double>>,
    BrownCase<Eigen::Vector<float, 2>, Eigen::Vector<float, 3>>,
    BrownCase<Eigen::VectorXf, Eigen::VectorXf>>;

OPTIMIST_TYPED_TEST_SUITE(Overdetermined, BrownTypes);

TYPED_TEST(Overdetermined, BadlyScaledZeroResidual) {
  using In     = typename TypeParam::In;
  using Out    = typename TypeParam::Out;
  using Scalar = typename TypeTrait<In>::Scalar;

  // The solution is (1e6, 2e-6): check the relative error componentwise
  Brown<In, Out> fun;
  Optimist::RootFinder::LevenbergMarquardt<In, Out> sol;
  sol.max_iterations(500);
  sol.max_function_evaluations(2000);
  sol.tolerance(solver_tolerance<Scalar>());
  In x_out;
  sol.rootfind(fun, fun.guess(0), x_out);
  EXPECT_TRUE(sol.converged());
  const Eigen::Vector<Scalar, Eigen::Dynamic> x(x_out);
  EXPECT_NEAR(x(0) / Scalar(1.0e6), 1, default_tolerance<Scalar>());
  EXPECT_NEAR(x(1) / Scalar(2.0e-6), 1, default_tolerance<Scalar>());
}

namespace {

  // Least-squares fit of y = c_0 exp(c_1 t) to eight perturbed samples of
  // 2 exp(-t/2): the residual at the solution is not zero
  class ExponentialFit
      : public DenseSystemFunction<Eigen::Vector2d,
                                   2,
                                   ExponentialFit,
                                   Eigen::Vector<double, 8>,
                                   8> {
   public:
    ExponentialFit() {
      for (Integer i{0}; i < 8; ++i) {
        m_t(i) = 0.5 * i;
        m_y(i) = 2 * std::exp(-0.5 * m_t(i)) + (i % 2 == 0 ? 0.05 : -0.05);
      }
      this->add_guess(this->make({1.0, 0.0}));
    }
    std::string name_impl() const {
      return "ExponentialFit";
    }
    void residual(const Dense &c, Dense &F) const {
      F = (c(0) * (c(1) * m_t.array()).exp() - m_y.array()).matrix();
    }
    void jacobian(const Dense &c, DenseMatrix &J) const {
      const Eigen::ArrayXd e{(c(1) * m_t.array()).exp()};
      J.col(0) = e.matrix();
      J.col(1) = (c(0) * m_t.array() * e).matrix();
    }
    void hessian(const Dense &c, const Integer i, DenseMatrix &H) const {
      const double e{std::exp(c(1) * m_t(i))};
      H << 0, m_t(i) * e, m_t(i) * e, c(0) * m_t(i) * m_t(i) * e;
    }

   private:
    Eigen::VectorXd m_t{8}, m_y{8};
  };

}  // namespace

TEST(Overdetermined, NonzeroResidual) {
  ExponentialFit fun;
  Optimist::RootFinder::LevenbergMarquardt<Eigen::Vector2d,
                                           Eigen::Vector<double, 8>>
      sol;
  sol.tolerance(1.0e-10);
  Eigen::Vector2d c;
  ASSERT_TRUE(sol.rootfind(fun, fun.guess(0), c));

  // Stationary point of |F|^2/2 with a nonzero residual
  Eigen::Vector<double, 8> F;
  Eigen::Matrix<double, 8, 2> J;
  ASSERT_TRUE(fun.evaluate(c, F));
  ASSERT_TRUE(fun.first_derivative(c, J));
  EXPECT_GT(F.norm(), 0.05);
  EXPECT_LT((J.transpose() * F).norm(), 1.0e-8);
  EXPECT_NEAR(c(0), 2.0, 0.1);
  EXPECT_NEAR(c(1), -0.5, 0.05);

  // Same minimizer as a general-purpose optimizer on |F|^2/2
  Optimist::Optimizer::BFGS<Eigen::Vector2d> bfgs;
  bfgs.tolerance(1.0e-9);
  Eigen::Vector2d c_bfgs;
  ASSERT_TRUE(bfgs.rootfind(fun, fun.guess(0), c_bfgs));
  EXPECT_LT((c - c_bfgs).norm(), 1.0e-6);
}
