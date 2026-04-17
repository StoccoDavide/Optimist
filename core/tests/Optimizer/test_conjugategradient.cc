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
    testing::Types<  // EllipticParaboloid<Eigen::Vector<float, 2>>,
                     // EllipticParaboloid<Eigen::Vector<float,
                     // Eigen::Dynamic>>,
                     // EllipticParaboloid<Eigen::SparseVector<float>>,
        Schaffer2<Eigen::Vector<double, 2>>,
        Schaffer2<Eigen::Vector<double, Eigen::Dynamic>>,
        Schaffer2<Eigen::SparseVector<double>>>;
// Schaffer2<double>,
// Brown<double>

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
  sol.verbose_mode(true);
  sol.tolerance(std::sqrt(Function::EPSILON));

  // Define methods to test
  using Method =
      typename Optimist::Optimizer::ConjugateGradient<Vector>::Method;
  auto methods = {
    // Method::FLETCHER_REEVES,
    Method::POLAK_RIBIERE,
    Method::POLAK_RIBIERE_PLUS
    // Method::HESTENES_STIEFEL,
    // Method::CONJUGATE_DESCENT,
    // Method::LIU_STOREY,
    // Method::DAI_YUAN,
    // Method::HAGER_ZHANG,
    // Method::HAGER_ZHANG_PLUS
  };

  // Test to change methods
  for (const auto &method : methods) {
    sol.method(method);

    Vector x_out;
    for (size_t i{0}; i < fun.guesses().size(); ++i) {
      sol.optimize(fun, fun.guess(i), x_out);
      EXPECT_TRUE(sol.converged());
      EXPECT_TRUE(fun.is_solution(x_out, std::cbrt(Function::EPSILON)));
    }
  }
}

//  sol.max_iterations(250);
//  sol.max_function_evaluations(5000);
//  sol.max_gradient_evaluations(5000);
//  sol.max_line_search_iterations(40);
//
