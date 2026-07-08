// /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *\
//  * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi. *
//  * *
//  * The Optimist project is distributed under the BSD 2-Clause License. *
//  * *
//  * Davide Stocco                                           Enrico Bertolazzi
//  *
//  * University of Trento                                 University of Trento
//  *
//  * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it
//  *
// \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// */

// // Optimist library
// #include "Optimist/Optimizer/NelderMead.hh"

// // Catch2 library
// #include <catch2/catch_template_test_macros.hpp>
// #include <catch2/catch_test_macros.hpp>
// #include <catch2/generators/catch_generators_range.hpp>

// #include "test_cost_functions.hh"

// using namespace Optimist;
// using namespace Optimist::TestSet;

// TEMPLATE_TEST_CASE("NelderMead", "[template]", TEST_COST_FUNCTIONS) {
//   TestType fun;
//   SECTION(fun.name()) {
//     Optimizer::NelderMead<double, fun.input_dimension()> sol;
//     SECTION(sol.name()) {
//       sol.task(fun.name());
//       sol.enable_verbose_mode();
//       typename TestType::InputType x_ini, x_out;
//       for (size_t i{0}; i < fun.guesses().size(); ++i) {
//         x_ini = fun.guess(i);
//         // Solve without damping
//         sol.disable_damped_mode();
//         sol.rootfind(fun, x_ini, x_out);
//         REQUIRE(sol.converged());
//         REQUIRE(fun.is_solution(x_out, TestType::SQRT_EPSILON));
//         // Solve with damping
//         sol.enable_damped_mode();
//         sol.rootfind(fun, x_ini, x_out);
//         REQUIRE(sol.converged());
//         REQUIRE(fun.is_solution(x_out, TestType::SQRT_EPSILON));
//       }
//     }
//   }
// }

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
#include "Optimist/Optimizer/NelderMead.hh"
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
                   // EllipticParaboloid<Eigen::SparseVector<float>>,
                   EllipticParaboloid<Eigen::Vector<double, 2>>,
                   EllipticParaboloid<Eigen::Vector<double, Eigen::Dynamic>>,
                   // EllipticParaboloid<Eigen::SparseVector<double>>,
                   Schaffer2<Eigen::Vector<float, 2>>,
                   Schaffer2<Eigen::Vector<float, Eigen::Dynamic>>,
                   // Schaffer2<Eigen::SparseVector<float>>,
                   Schaffer2<Eigen::Vector<double, 2>>,
                   Schaffer2<Eigen::Vector<double, Eigen::Dynamic>>
                   // Schaffer2<Eigen::SparseVector<double>>
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
  Optimist::Optimizer::NelderMead<Vector> sol;
  sol.task(fun.name());
  sol.verbose_mode(false);
  sol.tolerance(std::sqrt(Function::EPSILON));

  // Test to change methods
  Vector x_out;
  for (size_t i{0}; i < fun.guesses().size(); ++i) {
    sol.optimize(fun, fun.guess(i), x_out);
    EXPECT_TRUE(sol.converged());
    EXPECT_TRUE(fun.is_solution(x_out, 1.0e3 * std::cbrt(Function::EPSILON)));
  }
}
