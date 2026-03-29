/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#include "Optimist/TestSet/Brown.hh"
#include "Optimist/TestSet/EllipticParaboloid.hh"
#include "Optimist/TestSet/Schaffer2.hh"

#ifndef TEST_COST_FUNCTIONS
#define TEST_COST_FUNCTIONS                      \
  Optimist::TestSet::EllipticParaboloid<double>, \
      Optimist::TestSet::Schaffer2<double>, Optimist::TestSet::Brown<double>
#endif
