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

#ifndef INCLUDE_OPTIMIST_TESTSET_HH
#define INCLUDE_OPTIMIST_TESTSET_HH

#include "Optimist/TestSet/Booth.hh"
#include "Optimist/TestSet/Brown.hh"
#include "Optimist/TestSet/Cos.hh"
#include "Optimist/TestSet/Cosh.hh"
#include "Optimist/TestSet/EllipticParaboloid.hh"
#include "Optimist/TestSet/Linear.hh"
#include "Optimist/TestSet/Quadratic.hh"
#include "Optimist/TestSet/Rosenbrock.hh"
#include "Optimist/TestSet/Schaffer2.hh"
#include "Optimist/TestSet/Sin.hh"
#include "Optimist/TestSet/Sinh.hh"

namespace Optimist {

  /**
   * \brief Namespace for the Optimist library test set functions.
   */
  namespace TestSet {

    /**
     * List of cost functions.
     */
    static const std::vector<std::string> COST_FUNCTIONS = {
      "EllipticParaboloid",
      "Schaffer2"};

    /**
     * List of scalar-valued functions.
     */
    static const std::vector<std::string> SCALAR_FUNCTIONS =
        {"Linear", "Quadratic", "Sin", "Cos", "Sinh", "Cosh"};

    /**
     * List of vector-valued functions.
     */
    static const std::vector<std::string> VECTOR_FUNCTIONS = {
      "Booth",
      "Brown",
      "Rosenbrock",
    };

  }  // namespace TestSet

}  // namespace Optimist

#endif  // INCLUDE_OPTIMIST_TESTSET_HH