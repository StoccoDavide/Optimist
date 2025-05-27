/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimist library
#include "Optimist/RootFinder/Halley.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include "test_scalar_functions.hh"

using namespace Optimist;
using namespace Optimist::TestSet;

TEMPLATE_TEST_CASE("Halley", "[template]", TEST_SCALAR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    RootFinder::Halley<double> sol;
    sol.task(fun.name());
    typename TestType::InputType x_ini, x_out;
    for (Integer i{0}; i < static_cast<Integer>(fun.guesses().size()); ++i) {
      x_ini = fun.guess(i);
      // Solve without damping
      sol.disable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, TestType::EPSILON_LOW));
      // Solve with damping
      sol.enable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, TestType::EPSILON_LOW));
    }
  }
}
