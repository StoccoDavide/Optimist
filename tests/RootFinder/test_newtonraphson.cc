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
#include "Optimist/RootFinder/NewtonRaphson.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include "test_scalar_functions.hh"

using namespace Optimist;
using namespace Optimist::TestSet;

TEMPLATE_TEST_CASE("NewtonRaphson", "[template]", TEST_SCALAR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    RootFinder::NewtonRaphson<double> sol;
    sol.task(fun.name());
    typename TestType::InputType x_ini, x_out;
    for (size_t i{0}; i < fun.guesses().size(); ++i) {
      x_ini = fun.guess(i);
      // Solve without damping
      sol.disable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out,TestType::EPSILON_LOW));
      // Solve with damping
      sol.enable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out,TestType::EPSILON_LOW));
    }
  }
}
