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
#include "Optimist/RootFinder/Chandrupatla.hh"
#include "Optimist/TestSet.hh"

// Google Test framework
#include <gtest/gtest.h>

using namespace Optimist::TestSet;

// Type-parameterized test fixture for functions
template <typename FunctionType>
struct Functions : public testing::Test {
  using TestType = FunctionType;
};

using TestTypes = testing::Types<Linear<float>,
                                 Linear<double>,
                                 Quadratic<float>,
                                 Quadratic<double>,
                                 Cos<float>,
                                 Cos<double>,
                                 Sin<float>,
                                 Sin<double>,
                                 Sinh<float>,
                                 Sinh<double>>;

// Register the type-parameterized function tests
TYPED_TEST_SUITE(Functions, TestTypes);

// Test to solve functions
TYPED_TEST(Functions, Solve) {
  // Retrieve the function and scalar types
  using Function = TypeParam;
  using Scalar   = typename Function::Scalar;

  // Create function and solver instances
  Function fun;
  Optimist::RootFinder::Chandrupatla<Scalar> sol;
  sol.bounds(0.0, M_PI);
  sol.task(fun.name());
  sol.verbose_mode(true);
  sol.tolerance(std::sqrt(Function::EPSILON));

  Scalar x_out;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    sol.rootfind(fun, fun.guess(i), x_out);
    EXPECT_TRUE(sol.converged());
    EXPECT_TRUE(fun.is_solution(x_out, std::cbrt(Function::EPSILON)));
  }
}