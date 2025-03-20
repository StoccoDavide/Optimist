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

#include "Optimist/TestSet/CostFunction/EllipticParaboloid.hh"
#include "Optimist/TestSet/CostFunction/Schaffer2.hh"
#include "Optimist/TestSet/VectorFunction/Brown.hh"

#ifndef TEST_COST_FUNCTIONS
#define TEST_COST_FUNCTIONS                      \
  Optimist::TestSet::EllipticParaboloid<double>, \
  Optimist::TestSet::Schaffer2<double>,          \
  Optimist::TestSet::Brown<double>
#endif
