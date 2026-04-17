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
#include "Optimist/Optimizer/ConjugateGradient.hh"
#include "Optimist/TestSet.hh"

// Google Test framework
#include <gtest/gtest.h>

using namespace Optimist::TestSet;

// Type-parameterized test fixture for functions
template <typename FunctionType>
struct Functions : public testing::Test {
  using TestType = FunctionType;
};

using TestTypes =
    testing::Types<EllipticParaboloid<Eigen::Vector<float, 2>>,
                   EllipticParaboloid<Eigen::Vector<float, Eigen::Dynamic>>,
                   EllipticParaboloid<Eigen::SparseVector<float>>,
                   EllipticParaboloid<Eigen::Vector<double, 2>>,
                   EllipticParaboloid<Eigen::Vector<double, Eigen::Dynamic>>,
                   EllipticParaboloid<Eigen::SparseVector<double>>,
                   Schaffer2<Eigen::Vector<float, 2>>,
                   Schaffer2<Eigen::Vector<float, Eigen::Dynamic>>,
                   Schaffer2<Eigen::SparseVector<float>>,
                   Schaffer2<Eigen::Vector<double, 2>>,
                   Schaffer2<Eigen::Vector<double, Eigen::Dynamic>>,
                   Schaffer2<Eigen::SparseVector<double>>>;

// Register the type-parameterized function tests
TYPED_TEST_SUITE(Functions, TestTypes);

// Test to solve functions
TYPED_TEST(Functions, Solve) {
  // Retrieve the function and vector types
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  // Create function and solver instances
  Function fun;
  Optimist::Optimizer::ConjugateGradient<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(false);
  sol.tolerance(std::sqrt(Function::EPSILON));

  // Define methods to test
  using AlphaMethod =
      typename Optimist::Optimizer::ConjugateGradient<Vector>::AlphaMethod;
  using BetaMethod =
      typename Optimist::Optimizer::ConjugateGradient<Vector>::BetaMethod;
  auto alpha_methods = {AlphaMethod::ACCEPTED_STEP,
                        AlphaMethod::BARZILAI_BORWEIN};
  auto beta_methods  = {BetaMethod::FLETCHER_REEVES,
                        BetaMethod::POLAK_RIBIERE,
                        BetaMethod::POLAK_RIBIERE_PLUS,
                        BetaMethod::HESTENES_STIEFEL,
                        BetaMethod::CONJUGATE_DESCENT,
                        BetaMethod::LIU_STOREY,
                        BetaMethod::DAI_YUAN,
                        BetaMethod::HAGER_ZHANG,
                        BetaMethod::HAGER_ZHANG_PLUS};

  // Test to change methods
  for (const auto &alpha_method : alpha_methods) {
    for (const auto &beta_method : beta_methods) {
      sol.alpha_method(alpha_method);
      sol.beta_method(beta_method);

      Vector x_out;
      for (size_t i{0}; i < fun.guesses().size(); ++i) {
        sol.optimize(fun, fun.guess(i), x_out);
        EXPECT_TRUE(sol.converged());
        EXPECT_TRUE(
            fun.is_solution(x_out, 1.0e3 * std::cbrt(Function::EPSILON)));
      }
    }
  }
}
