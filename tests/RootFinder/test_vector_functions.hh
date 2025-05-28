/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#include "Optimist/TestSet/Booth.hh"
#include "Optimist/TestSet/Rosenbrock.hh"
#include "Optimist/TestSet/Test11.hh"

#ifndef TEST_VECTOR_FUNCTIONS
#define TEST_VECTOR_FUNCTIONS             \
  Optimist::TestSet::Test11<double>,      \
  Optimist::TestSet::Booth<double>,       \
  Optimist::TestSet::Rosenbrock2<double>, \
  Optimist::TestSet::Rosenbrock4<double>
#endif
