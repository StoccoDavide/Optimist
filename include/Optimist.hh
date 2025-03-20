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

#ifndef OPTIMIST_DEFAULT_INTEGER_TYPE
#define OPTIMIST_DEFAULT_INTEGER_TYPE int
#endif

/**
* \brief Namespace for the Optimist library.
*/
namespace Optimist
{

  /**
  * \brief The Integer type as used for the API.
  *
  * The Integer type, \c \#define the preprocessor symbol \c OPTIMIST_DEFAULT_INTEGER_TYPE. The default
  * value is \c int.
  */
  using Integer = OPTIMIST_DEFAULT_INTEGER_TYPE;

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
    for (Integer i{0}; i < N; ++i) {line += table_horizontal_line();}
    return line;
  }

  /**
  * \brief Retrieve the Unicode character for the vertical line of a table.
  * \return Unicode character for the vertical line of a table.
  */
  static inline std::string table_vertical_line() {return "│";}

} // namespace Optimist

#endif // INCLUDE_OPTIMIST_HH
