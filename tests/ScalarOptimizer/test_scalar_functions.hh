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

#include "Optimist/TestSet/ScalarFunction/Quadratic.hh"
#include "Optimist/TestSet/ScalarFunction/Cos.hh"
#include "Optimist/TestSet/ScalarFunction/Sin.hh"
#include "Optimist/TestSet/ScalarFunction/Cosh.hh"

#ifndef TEST_SCALAR_FUNCTIONS
#define TEST_SCALAR_FUNCTIONS           \
  Optimist::TestSet::Quadratic<double>, \
  Optimist::TestSet::Cos<double>,       \
  Optimist::TestSet::Sin<double>,       \
  Optimist::TestSet::Cosh<double>
#endif
