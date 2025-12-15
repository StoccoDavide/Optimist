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
#include "Optimist/RootFinder/Greenstadt.hh"

// Google Test framework
#include <gtest/gtest.h>

using namespace Optimist::TestSet;

// Type-parameterized test fixture for functions
template <typename FunctionType>
struct Functions : public testing::Test {
    using TestType = FunctionType;
};

using TestTypes = testing::Types<
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

// Register the type-parameterized function tests
TYPED_TEST_SUITE(Functions, TestTypes);

// Test to solve functions
TYPED_TEST(Functions, Solve) {

  // Retrieve the function and vector types
  using Function = TypeParam;
  using Vector   = typename Function::VectorTrait::Type;

  // Create function and solver instances
  Function fun;
  Optimist::RootFinder::Greenstadt<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(false);
  sol.tolerance(std::sqrt(Function::EPSILON));

  // Define methods to test
  using Method = typename Optimist::RootFinder::Greenstadt<Vector>::Method;
  auto methods = {Method::ONE, Method::TWO};

  // Test to change methods
  for (const auto & method : methods) {
    sol.method(method);

    Vector x_out;
    for (size_t i{0}; i < fun.guesses().size(); ++i) {

      // Solve without damping
      sol.disable_damped_mode();
      sol.rootfind(fun, fun.guess(i), x_out);
      EXPECT_TRUE(sol.converged());
      EXPECT_TRUE(fun.is_solution(x_out, std::cbrt(Function::EPSILON)));

      // Solve with damping
      sol.enable_damped_mode();
      sol.rootfind(fun, fun.guess(i), x_out);
      EXPECT_TRUE(sol.converged());
      EXPECT_TRUE(fun.is_solution(x_out, std::cbrt(Function::EPSILON)));
    }
  }
}
