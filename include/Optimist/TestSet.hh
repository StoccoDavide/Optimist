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

#ifndef INCLUDE_OPTIMIST_TESTSET_HH
#define INCLUDE_OPTIMIST_TESTSET_HH

// Optimist library
#include "Optimist.hh"

namespace Optimist
{
  /*\
   |   _____         _   ____       _
   |  |_   _|__  ___| |_/ ___|  ___| |_
   |    | |/ _ \/ __| __\___ \ / _ \ __|
   |    | |  __/\__ \ |_ ___) |  __/ |_
   |    |_|\___||___/\__|____/ \___|\__|
   |
  \*/

  /**
  * \brief Namespace for the Optimist library test set functions.
  */
  namespace TestSet
  {
    static constexpr Real PI{3.14159265358979323846}; /**< The value of \f$\pi\f$. */

    /**
    * Map of vector-valued functions.
    */
    static const std::vector<std::string> vector_functions = {
      "Wood"
    };

  } // namespace TestSet

} // namespace Optimist

// Objective functions
#include "TestSet/CostFunction/Schaffer2.hxx"

// Scalar-valued functions
#include "TestSet/ScalarFunction/Sin.hxx"
#include "TestSet/ScalarFunction/Cos.hxx"
#include "TestSet/ScalarFunction/Cosh.hxx"

// Vector-valued functions
#include "TestSet/VectorFunction/Booth.hxx"
#include "TestSet/VectorFunction/Brown.hxx"
#include "TestSet/VectorFunction/Rosenbrock2.hxx"

#endif // INCLUDE_OPTIMIST_TESTSET_HH