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
#include "Optimist/TestSet/Sin.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh"
#include "Optimist/RootFinder/Greenstadt.hh"

// Run all the tests.
int main() {

  Optimist::TestSet::Sin<double> fun;
  Optimist::RootFinder::NewtonRaphson<double> sol;
  double x_out;
  sol.rootfind(fun, 0.0, x_out);
  Optimist::RootFinder::Greenstadt<double, 4> greenstadt;

  return 1;
}
