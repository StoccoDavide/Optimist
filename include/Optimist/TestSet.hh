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

    /**
    * Print Optimist library test-set information on a string.
    * \return A string with the Optimist library test-set information.
    */
    std::string TestSetInfo() {
      std::ostringstream os;
      os
        << Optimist::Info() << std::endl
        << " Test set cost functions:" << std::endl
        << std::accumulate(COST_FUNCTIONS.begin(), COST_FUNCTIONS.end(), std::string(" ")).substr(0, 100) << std::endl
        << " Test set scalar-valued functions:" << std::endl
        << std::accumulate(SCALAR_FUNCTIONS.begin(), SCALAR_FUNCTIONS.end(), std::string(" ")).substr(0, 100) << std::endl
        << " Test set vector-valued functions:" << std::endl
        << std::accumulate(VECTOR_FUNCTIONS.begin(), VECTOR_FUNCTIONS.end(), std::string(" ")).substr(0, 100) << std::endl;

      return os.str();
    }

    /**
    * Print Optimist library test-set information on a stream.
    * \param[in] os Output stream.
    */
    void TestSetInfo(std::ostream &os) {os << Info();}



  } // namespace TestSet

} // namespace Optimist

#endif // INCLUDE_OPTIMIST_TESTSET_HH