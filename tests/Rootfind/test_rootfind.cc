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
#include "Optimist.hh"
#include "Optimist/TestSet.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace Optimist;
using namespace Optimist::TestSet;

#define TEST_SCALAR_FUNCTIONS Cos, Cosh, Sin
#define TEST_VECTOR_FUNCTIONS Booth, Rosenbrock<2>, Rosenbrock<4>, Rosenbrock<8>

TEMPLATE_TEST_CASE("Scalar - Newton", "[template]", TEST_SCALAR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    ScalarRootFinder::Newton sol;
    sol.task(fun.name());
    typename TestType::InputType x_ini, x_out;
    for (Integer i{0}; i < static_cast<Integer>(fun.guesses().size()); ++i) {
      x_ini = fun.guess(i);
      // Solve without damping
      sol.disable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
      // Solve with damping
      sol.enable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
    }
  }
}

TEMPLATE_TEST_CASE("Vector - Newton", "[template]", TEST_VECTOR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    RootFinder::Newton<fun.input_dimension()> sol;
    sol.task(fun.name());
    typename TestType::InputType x_ini, x_out;
    for (Integer i{0}; i < static_cast<Integer>(fun.guesses().size()); ++i) {
      x_ini = fun.guess(i);
      // Solve without damping
      sol.disable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
      // Solve with damping
      sol.enable_damped_mode();
      sol.rootfind(fun, x_ini, x_out);
      REQUIRE(sol.converged());
      REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
    }
  }
}

TEMPLATE_TEST_CASE("Vector - Broyden", "[template]", TEST_VECTOR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    RootFinder::Broyden<fun.input_dimension()> sol;
    using BroydenMethod = typename RootFinder::Broyden<fun.input_dimension()>::Method;
    auto met = GENERATE(BroydenMethod::GOOD, BroydenMethod::BAD, BroydenMethod::COMBINED);
    sol.method(met);
    SECTION(sol.name()) {
      sol.task(fun.name());
      typename TestType::InputType x_ini, x_out;
      for (Integer i{0}; i < static_cast<Integer>(fun.guesses().size()); ++i) {
        x_ini = fun.guess(i);
        // Solve without damping
        sol.disable_damped_mode();
        sol.rootfind(fun, x_ini, x_out);
        REQUIRE(sol.converged());
        REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
        // Solve with damping
        sol.enable_damped_mode();
        sol.rootfind(fun, x_ini, x_out);
        REQUIRE(sol.converged());
        REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
      }
    }
  }
}

TEMPLATE_TEST_CASE("Vector - Greenstadt", "[template]", TEST_VECTOR_FUNCTIONS) {
  TestType fun;
  SECTION(fun.name()) {
    RootFinder::Greenstadt<fun.input_dimension()> sol;
    using GreenstadtMethod = typename RootFinder::Greenstadt<fun.input_dimension()>::Method;
    auto met = GENERATE(GreenstadtMethod::ONE, GreenstadtMethod::TWO);
    sol.method(met);
    SECTION(sol.name()) {
      sol.task(fun.name());
      typename TestType::InputType x_ini, x_out;
      for (Integer i{0}; i < static_cast<Integer>(fun.guesses().size()); ++i) {
        x_ini = fun.guess(i);
        // Solve without damping
        sol.disable_damped_mode();
        sol.rootfind(fun, x_ini, x_out);
        REQUIRE(sol.converged());
        REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
        // Solve with damping
        sol.enable_damped_mode();
        sol.rootfind(fun, x_ini, x_out);
        REQUIRE(sol.converged());
        REQUIRE(fun.is_solution(x_out, EPSILON_LOW));
      }
    }
  }
}
