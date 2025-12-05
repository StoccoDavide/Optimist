/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_SCALAR_FUNCTIONS_HH
#define OPTIMIST_TESTSET_SCALAR_FUNCTIONS_HH

#include "Optimist/TestSet/Quadratic.hh"
#include "Optimist/TestSet/Cos.hh"
#include "Optimist/TestSet/Sin.hh"
#include "Optimist/TestSet/Cosh.hh"

#ifndef TEST_SCALAR_FUNCTIONS
#define TEST_SCALAR_FUNCTIONS           \
  Optimist::TestSet::Quadratic<double>, \
  Optimist::TestSet::Cos<double>,       \
  Optimist::TestSet::Sin<double>,       \
  Optimist::TestSet::Cosh<double>
#endif

#endif /* OPTIMIST_TESTSET_SCALAR_FUNCTIONS_HH */
