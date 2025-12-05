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
//#include "Optimist/RootFinder/NelderMead.hh"    // For testing cost functions

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
  Quadratic<float>,  Cos<float>,  Sin<float>, Cosh<float>,
  Quadratic<double>, Cos<double>, Sin<double>, Cosh<double>
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
  sol.tolerance(1.0e3*Function::EPSILON_LOW); // Very loose, we just want to test function evaluation

  bool converged_at_least_once = false;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    Scalar x_ini{fun.guess(i)}, x_out;
    sol.rootfind(fun, x_ini, x_out);

    // Check convergence
    if (sol.converged() && fun.is_solution(x_out, Function::EPSILON_LOW)) {
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
  //Booth<Eigen::Vector<float, 2>>,
  //Booth<Eigen::Vector<float, Eigen::Dynamic>>,
  //Booth<Eigen::SparseVector<float>>,
  Booth<Eigen::Vector<double, 2>>,
  Booth<Eigen::Vector<double, Eigen::Dynamic>>
  //Booth<Eigen::SparseVector<double>>
>;

// Register the type-parameterized vector function tests
TYPED_TEST_SUITE(VectorFunctions, VectorTestTypes);

// Test to solve vector functions
TYPED_TEST(VectorFunctions, Solve) {
  // Retrieve the function and vector types
  using Function = TypeParam;
  using Vector   = typename Optimist::RetriveType<Function>::First;

  // Create function and solver instances
  Function fun;
  Optimist::RootFinder::Newton<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(true);
  sol.tolerance(1.0e3*Function::EPSILON_LOW); // Very loose, we just want to test function evaluation

  bool converged_at_least_once = false;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    Vector x_ini{fun.guess(i)}, x_out;
    sol.rootfind(fun, x_ini, x_out);

    // Check convergence
    if (sol.converged() && fun.is_solution(x_out, Function::EPSILON_LOW)) {
      converged_at_least_once = true;
    }
  }
  ASSERT_TRUE(converged_at_least_once);
}