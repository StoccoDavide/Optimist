/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist.hh"
#include "Optimist/TestSet.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh" // For testing scalar functions
#include "Optimist/RootFinder/Newton.hh"        // For testing vector functions
#include "Optimist/Optimizer/NelderMead.hh"    // For testing cost functions

// Test functions for the Optimist library
// Google Test framework
#include <gtest/gtest.h>

using namespace Optimist::TestSet;

// Type-parameterized test fixture for scalar functions
template <typename FunctionType>
struct ScalarFunctions : public testing::Test {
    using TestType = FunctionType;
};

using ScalarTestTypes = testing::Types<
  Linear<float>,
  Linear<double>,
  Quadratic<float>,
  Quadratic<double>,
  Cos<float>,
  Cos<double>,
  Sin<float>,
  Sin<double>,
  Sinh<float>,
  Sinh<double>
>;

// Register the type-parameterized scalar function tests
TYPED_TEST_SUITE(ScalarFunctions, ScalarTestTypes);

// Test to solve scalar functions
TYPED_TEST(ScalarFunctions, Solve) {

  // Retrieve the function and scalar types
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  // Create function and solver instances
  Function fun;
  Optimist::RootFinder::NewtonRaphson<Scalar> sol;
  sol.task(fun.name());
  sol.verbose_mode(false);
  sol.tolerance(std::sqrt(Function::EPSILON));

  bool converged_at_least_once{false};
  Scalar x_out;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    sol.rootfind(fun, fun.guess(i), x_out);

    // Check convergence
    if (sol.converged() && fun.is_solution(x_out, std::cbrt(Function::EPSILON))) {
      converged_at_least_once = true;
    }
  }
  ASSERT_TRUE(converged_at_least_once);
}

// Type-parameterized test fixture for vector functions
template <typename FunctionType>
struct VectorFunctions : public testing::Test {
    using TestType = FunctionType;
};

using VectorTestTypes = testing::Types<
  Linear1<Eigen::Vector<float, 1>>,
  Linear1<Eigen::Vector<float, Eigen::Dynamic>>,
  Linear1<Eigen::SparseVector<float>>,
  Linear1<Eigen::Vector<double, 1>>,
  Linear1<Eigen::Vector<double, Eigen::Dynamic>>,
  Linear1<Eigen::SparseVector<double>>,
  Booth<Eigen::Vector<float, 2>>,
  Booth<Eigen::Vector<float, Eigen::Dynamic>>,
  Booth<Eigen::SparseVector<float>>,
  Booth<Eigen::Vector<double, 2>>,
  Booth<Eigen::Vector<double, Eigen::Dynamic>>,
  Booth<Eigen::SparseVector<double>>,
  Rosenbrock2<Eigen::Vector<float, 2>>,
  Rosenbrock2<Eigen::Vector<float, Eigen::Dynamic>>,
  Rosenbrock2<Eigen::SparseVector<float>>,
  Rosenbrock2<Eigen::Vector<double, 2>>,
  Rosenbrock2<Eigen::Vector<double, Eigen::Dynamic>>,
  Rosenbrock2<Eigen::SparseVector<double>>,
  Rosenbrock4<Eigen::Vector<float, 4>>,
  Rosenbrock4<Eigen::Vector<float, Eigen::Dynamic>>,
  Rosenbrock4<Eigen::SparseVector<float>>,
  Rosenbrock4<Eigen::Vector<double, 4>>,
  Rosenbrock4<Eigen::Vector<double, Eigen::Dynamic>>,
  Rosenbrock4<Eigen::SparseVector<double>>,
  Rosenbrock6<Eigen::Vector<float, 6>>,
  Rosenbrock6<Eigen::Vector<float, Eigen::Dynamic>>,
  Rosenbrock6<Eigen::SparseVector<float>>,
  Rosenbrock6<Eigen::Vector<double, 6>>,
  Rosenbrock6<Eigen::Vector<double, Eigen::Dynamic>>,
  Rosenbrock6<Eigen::SparseVector<double>>,
  Rosenbrock8<Eigen::Vector<float, 8>>,
  Rosenbrock8<Eigen::Vector<float, Eigen::Dynamic>>,
  Rosenbrock8<Eigen::SparseVector<float>>,
  Rosenbrock8<Eigen::Vector<double, 8>>,
  Rosenbrock8<Eigen::Vector<double, Eigen::Dynamic>>,
  Rosenbrock8<Eigen::SparseVector<double>>
>;

// Register the type-parameterized vector function tests
TYPED_TEST_SUITE(VectorFunctions, VectorTestTypes);

// Test to solve vector functions
TYPED_TEST(VectorFunctions, Solve) {
  // Retrieve the function and vector types
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  // Create function and solver instances
  Function fun;
  Optimist::RootFinder::Newton<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(false);
  sol.tolerance(std::sqrt(Function::EPSILON));

  bool converged_at_least_once{false};
  Vector x_out;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    sol.rootfind(fun, fun.guess(i), x_out);

    // Check convergence
    if (sol.converged() && fun.is_solution(x_out, std::cbrt(Function::EPSILON))) {
      converged_at_least_once = true;
    }
  }
  ASSERT_TRUE(converged_at_least_once);
}

// Type-parameterized test fixture for cost functions
template <typename FunctionType>
struct CostFunctions : public testing::Test {
    using TestType = FunctionType;
};

// Type-parameterized test fixture for cost functions
using CostTestTypes = testing::Types<
  //Brown<Eigen::Vector<float, 2>, Eigen::Vector<float, 3>>,
  //Brown<Eigen::Vector<float, Eigen::Dynamic>, Eigen::Vector<float, Eigen::Dynamic>>,
  //Brown<Eigen::SparseVector<float>, Eigen::SparseVector<float>>,
  //Brown<Eigen::Vector<double, 2>, Eigen::Vector<double, 3>>,
  //Brown<Eigen::Vector<double, Eigen::Dynamic>, Eigen::Vector<double, Eigen::Dynamic>>,
  //Brown<Eigen::SparseVector<double>, Eigen::SparseVector<double>>,
  Schaffer2<Eigen::Vector<float, 2>>
  //Schaffer2<Eigen::Vector<float, Eigen::Dynamic>>,
  //Schaffer2<Eigen::SparseVector<float>>,
  //Schaffer2<Eigen::Vector<double, 2>>,
  //Schaffer2<Eigen::Vector<double, Eigen::Dynamic>>,
  //Schaffer2<Eigen::SparseVector<double>>,
  //EllipticParaboloid<Eigen::Vector<float, 2>>,
  //EllipticParaboloid<Eigen::Vector<float, Eigen::Dynamic>>,
  //EllipticParaboloid<Eigen::SparseVector<float>>,
  //EllipticParaboloid<Eigen::Vector<double, 2>>,
  //EllipticParaboloid<Eigen::Vector<double, Eigen::Dynamic>>,
  //EllipticParaboloid<Eigen::SparseVector<double>>
>;

// Register the type-parameterized cost function tests
TYPED_TEST_SUITE(CostFunctions, CostTestTypes);

// Test to solve cost functions
TYPED_TEST(CostFunctions, SolveCost) {
  // Retrieve the function and vector types
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  // Create function and solver instances
  Function fun;
  Optimist::Optimizer::NelderMead<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(true);
  sol.tolerance(std::sqrt(Function::EPSILON));

  bool converged_at_least_once{false};
  Vector x_out;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    sol.optimize(fun, fun.guess(i), x_out);

    // Check convergence
    if (sol.converged() && fun.is_solution(x_out, std::cbrt(Function::EPSILON))) {
      converged_at_least_once = true;
    }
  }
  ASSERT_TRUE(converged_at_least_once);
}

