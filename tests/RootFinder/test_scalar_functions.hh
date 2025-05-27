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

#include "Optimist/TestSet/Linear.hh"
#include "Optimist/TestSet/Quadratic.hh"
#include "Optimist/TestSet/Cos.hh"
#include "Optimist/TestSet/Sin.hh"
#include "Optimist/TestSet/Sinh.hh"

#ifndef TEST_SCALAR_FUNCTIONS
#define TEST_SCALAR_FUNCTIONS           \
  Optimist::TestSet::Linear<double>,    \
  Optimist::TestSet::Quadratic<double>, \
  Optimist::TestSet::Cos<double>,       \
  Optimist::TestSet::Sin<double>,       \
  Optimist::TestSet::Sinh<double>
#endif
