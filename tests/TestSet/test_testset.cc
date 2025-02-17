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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators_range.hpp>

// Run all the tests.
int main() {

  Optimist::Info(std::cout);
  Optimist::TestSet::TestSetInfo(std::cout);
  Optimist::TestSet::Cos cos;
  Optimist::ScalarRootFinder::Newton sol;
  double x_out;
  sol.rootfind(cos, 0.0, x_out);
  Optimist::RootFinder::Greenstadt<4> greenstadt;

  return 1;
}
