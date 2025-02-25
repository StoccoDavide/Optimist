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

#ifndef INCLUDE_OPTIMIST_HH
#define INCLUDE_OPTIMIST_HH

// C++ standard libraries
#include <iostream>
#include <ios>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <numeric>
#include <algorithm>

// Eigen library
#include <Eigen/Dense>

// Print Optimist errors
#ifndef OPTIMIST_ERROR
#define OPTIMIST_ERROR(MSG)             \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for Optimist
#ifndef OPTIMIST_ASSERT
#define OPTIMIST_ASSERT(COND, MSG) \
  if (!(COND))                     \
  {                                \
    OPTIMIST_ERROR(MSG);           \
  }
#endif

// Warning for Optimist
#ifndef OPTIMIST_WARNING
#define OPTIMIST_WARNING(MSG)       \
  {                                 \
    std::cout << MSG << std::endl;  \
  }
#endif

// Warning assert for Optimist
#ifndef OPTIMIST_ASSERT_WARNING
#define OPTIMIST_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                             \
  {                                        \
    OPTIMIST_WARNING(MSG);                 \
  }
#endif


// Define the basic constants for Optimist
#ifndef OPTIMIST_BASIC_CONSTANTS
#define OPTIMIST_BASIC_CONSTANTS(Real) \
  static constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};     /**< Machine epsilon epsilon static constant value. */ \
  static constexpr Real EPSILON_HIGH{1.0e-12};                             /**< High precision epsilon static constant value. */ \
  static constexpr Real EPSILON_MEDIUM{1.0e-10};                           /**< Medium precision epsilon static constant value. */ \
  static constexpr Real EPSILON_LOW{1.0e-08};                              /**< Low precision epsilon static constant value. */ \
  static constexpr Real INFTY{std::numeric_limits<Real>::infinity()};      /**< Infinity static constant value. */ \
  static constexpr Real QUIET_NAN{std::numeric_limits<Real>::quiet_NaN()}; /**< Not-a-number static constant value. */
#endif

/**
* \brief Namespace for the Optimist library.
*/
namespace Optimist
{

  using Integer = int; /**< Integer number type. */

  /**
  * \brief Retrieve the Unicode character for the top-left corner of a table.
  * \return Unicode character for the top-left corner of a table.
  */
  static inline std::string table_top_left_corner() {return "┌";}

  /**
  * \brief Retrieve the Unicode character for the top-right corner of a table.
  * \return Unicode character for the top-right corner of a table.
  */
  static inline std::string table_top_right_corner() {return "┐";}

  /**
  * \brief Retrieve the Unicode character for the bottom-left corner of a table.
  * \return Unicode character for the bottom-left corner of a table.
  */
  static inline std::string table_bottom_left_corner() {return "└";}

  /**
  * \brief Retrieve the Unicode character for the bottom-right corner of a table.
  * \return Unicode character for the bottom-right corner of a table.
  */
  static inline std::string table_bottom_right_corner() {return "┘";}

  /**
  * \brief Retrieve the Unicode character for the left junction of a table.
  * \return Unicode character for the left junction of a table.
  */
  static inline std::string table_left_junction() {return "├";}

  /**
  * \brief Retrieve the Unicode character for the right junction of a table.
  * \return Unicode character for the right junction of a table.
  */
  static inline std::string table_right_junction() {return "┤";}

  /**
  * \brief Retrieve the Unicode character for the top junction of a table.
  * \return Unicode character for the top junction of a table.
  */
  static inline std::string table_top_junction() {return "┬";}

  /**
  * \brief Retrieve the Unicode character for the bottom junction of a table.
  * \return Unicode character for the bottom junction of a table.
  */
  static inline std::string table_bottom_junction() {return "┴";}

  /**
  * \brief Retrieve the Unicode character for the center cross of a table.
  * \return Unicode character for the center cross of a table.
  */
  static inline std::string table_center_cross() {return "┼";}

  /**
  * \brief Retrieve the Unicode character for the horizontal line of a table.
  * \return Unicode character for the horizontal line of a table.
  */
  static inline std::string table_horizontal_line() {return "─";}

  /**
  * \brief Retrieve the Unicode character for a number of horizontal lines of a table.
  * \return Unicode character for the horizontal lines of a table.
  * \tparam N Number of horizontal lines.
  */
  template <Integer N>
  static inline std::string table_horizontal_line() {
    std::string line;
    for (Integer i = 0; i < N; ++i) {line += table_horizontal_line();}
    return line;
  }

  /**
  * \brief Retrieve the Unicode character for the vertical line of a table.
  * \return Unicode character for the vertical line of a table.
  */
  static inline std::string table_vertical_line() {return "│";}

  /**
  * Print Optimist library information on a string.
  * \return A string with the Optimist library information.
  */
  inline std::string Info() {
    std::ostringstream os;
    os
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl
      << "* Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *" << std::endl
      << "*                                                                                               *" << std::endl
      << "* The Optimist project is distributed under the BSD 2-Clause License.                           *" << std::endl
      << "*                                                                                               *" << std::endl
      << "* Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *" << std::endl
      << "* University of Trento               University of Trento                  University of Trento *" << std::endl
      << "* davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *" << std::endl
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    return os.str();
  }

  /**
  * Print Optimist library information on a stream.
  * \param[in] os Output stream.
  */
  inline void Info(std::ostream &os) {os << Info();}

} // namespace Optimist

// The functions
#include "Optimist/Function.hxx"
#include "Optimist/Function/ScalarFunction.hxx"
#include "Optimist/Function/VectorFunction.hxx"
#include "Optimist/Function/CostFunction.hxx"

// The solver base class
#include "Optimist/Solver.hxx"

// Scalar root-finding solvers
#include "Optimist/ScalarRootFinder.hxx"
#include "Optimist/ScalarRootFinder/Bracketing.hxx"
#include "Optimist/ScalarRootFinder/Algo748.hxx"
#include "Optimist/ScalarRootFinder/Chandrupatla.hxx"
#include "Optimist/ScalarRootFinder/Chebyshev.hxx"
#include "Optimist/ScalarRootFinder/Halley.hxx"
#include "Optimist/ScalarRootFinder/Newton.hxx"
#include "Optimist/ScalarRootFinder/Varona.hxx"

// Root-finding solvers
#include "Optimist/RootFinder.hxx"
#include "Optimist/RootFinder/Newton.hxx"
#include "Optimist/RootFinder/QuasiNewton.hxx"
#include "Optimist/RootFinder/Broyden.hxx"
#include "Optimist/RootFinder/Greenstadt.hxx"

// Optimization solvers
#include "Optimist/Optimizer.hxx"
#include "Optimist/Optimizer/NelderMead.hxx"

// Scalar optimization solvers
#include "Optimist/ScalarOptimizer.hxx"

#endif // INCLUDE_OPTIMIST_HH
