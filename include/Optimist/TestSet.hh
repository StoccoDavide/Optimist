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

// Cost functions
#include "TestSet/CostFunction/EllipticParaboloid.hxx"
#include "TestSet/CostFunction/Schaffer2.hxx"

// Scalar-valued functions
#include "TestSet/ScalarFunction/Linear.hxx"
#include "TestSet/ScalarFunction/Quadratic.hxx"
#include "TestSet/ScalarFunction/Sin.hxx"
#include "TestSet/ScalarFunction/Cos.hxx"
#include "TestSet/ScalarFunction/Sinh.hxx"
#include "TestSet/ScalarFunction/Cosh.hxx"

// Vector-valued functions
#include "TestSet/VectorFunction/Booth.hxx"
#include "TestSet/VectorFunction/Brown.hxx"
#include "TestSet/VectorFunction/Rosenbrock.hxx"

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

    /**
    * List of cost functions.
    */
    static const std::vector<std::string> COST_FUNCTIONS = {
      "Schaffer2"
    };

    /**
    * List of scalar-valued functions.
    */
    static const std::vector<std::string> SCALAR_FUNCTIONS = {
      "Sin",
      "Cos",
      "Cosh"
    };

    /**
    * List of vector-valued functions.
    */
    static const std::vector<std::string> VECTOR_FUNCTIONS = {
      "Booth",
      "Brown",
      "Rosenbrock",
    };

  } // namespace TestSet

} // namespace Optimist

#endif // INCLUDE_OPTIMIST_TESTSET_HH