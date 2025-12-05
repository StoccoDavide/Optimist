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

#ifndef INCLUDE_OPTIMIST_HH
#define INCLUDE_OPTIMIST_HH

// C++ standard libraries
#include <iostream>
#include <ios>
#include <iomanip>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
#include <map>
#include <memory>
#include <numeric>
#include <algorithm>
#include <type_traits>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

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
#define OPTIMIST_BASIC_CONSTANTS(Scalar) \
  static constexpr Scalar EPSILON{std::numeric_limits<Scalar>::epsilon()};     /**< Machine epsilon epsilon static constant value. */ \
  static constexpr Scalar EPSILON_HIGH{1.0e-12};                             /**< High precision epsilon static constant value. */ \
  static constexpr Scalar EPSILON_MEDIUM{1.0e-10};                           /**< Medium precision epsilon static constant value. */ \
  static constexpr Scalar EPSILON_LOW{1.0e-08};                              /**< Low precision epsilon static constant value. */ \
  static constexpr Scalar INFTY{std::numeric_limits<Scalar>::infinity()};      /**< Infinity static constant value. */ \
  static constexpr Scalar QUIET_NAN{std::numeric_limits<Scalar>::quiet_NaN()}; /**< Not-a-number static constant value. */
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

  /*\
   |   _____                _____          _ _
   |  |_   _|   _ _ __   __|_   _| __ __ _(_) |_ ___
   |    | || | | | '_ \ / _ \| || '__/ _` | | __/ __|
   |    | || |_| | |_) |  __/| || | | (_| | | |_\__ \
   |    |_| \__, | .__/ \___||_||_|  \__,_|_|\__|___/
   |        |___/|_|
  \*/

  /**
   * Traits class for vectors and matrices (fallback for unsupported types).
   * \tparam T The type to be specialized.
   */
  template <typename T, typename Enable = void>
  struct TypeTrait;

  /**
   * Traits class for scalar types.
   * \tparam Scalar The scalar type.
   */
  template <typename ScalarType>
  struct TypeTrait<ScalarType, std::enable_if_t<std::is_floating_point<ScalarType>::value>>
  {
    using Scalar = ScalarType;
    static constexpr Integer Dimension{1};
    static constexpr bool IsScalar{true};
    static constexpr bool IsEigen{false};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{false};
  };

  /**
   * Traits class for fixed-size dense Eigen column vectors.
   * \tparam Scalar The scalar type.
   * \tparam N The vector dimension.
   */
  template <typename ScalarType, Integer N>
  struct TypeTrait<Eigen::Matrix<ScalarType, N, 1>, std::enable_if_t<(N > 0)>>
  {
    using Scalar = ScalarType;
    static constexpr Integer Dimension{N};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{true};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{false};
  };

  /**
   * Traits class for dynamic-size dense Eigen column vectors.
   * \tparam Scalar The scalar type.
   */
  template <typename ScalarType>
  struct TypeTrait<Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>>
  {
    using Scalar = ScalarType;
    static constexpr Integer Dimension{Eigen::Dynamic};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{true};
    static constexpr bool IsSparse{false};
  };

  /**
   * Traits class for sparse Eigen column vectors.
   * \tparam Scalar The scalar type.
   */
  template <typename ScalarType>
  struct TypeTrait<Eigen::SparseVector<ScalarType>>
  {
    using Scalar = ScalarType;
    static constexpr Integer Dimension{Eigen::Dynamic};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{true};
  };

  /**
   * Traits class for sparse Eigen matrices.
   * \tparam T The sparse Eigen matrix type.
   */
  template <typename ScalarType, Integer Options, typename StorageIndex>
  struct TypeTrait<Eigen::SparseMatrix<ScalarType, Options, StorageIndex>>
  {
    using Scalar = ScalarType;
    static constexpr Integer Dimension{Eigen::Dynamic};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{true};
  };

  /**
   * Traits class for understanding if two types are the comaptible (fallback).
   * \tparam FirstType The first type.
   * \tparam SecondType The second type.
   */
  template <typename FirstType, typename SecondType>
  struct AreCompatibleTypes;

  /**
   * Traits class for understanding if two types are the comaptible.
   * Specialization for compatible types.
   * \tparam FirstType The first type.
   * \tparam SecondType The second type.
   */
  template <typename FirstType, typename SecondType>
  struct AreCompatibleTypes<TypeTrait<FirstType>, TypeTrait<SecondType>>
  {
    static constexpr bool Value{
      (TypeTrait<FirstType>::IsScalar && TypeTrait<SecondType>::IsScalar) ||
      ((TypeTrait<FirstType>::IsEigen && TypeTrait<SecondType>::IsEigen) &&
      ((TypeTrait<FirstType>::IsFixed && TypeTrait<SecondType>::IsFixed &&
        (TypeTrait<FirstType>::Dimension == TypeTrait<SecondType>::Dimension)) ||
      (TypeTrait<FirstType>::IsDynamic && TypeTrait<SecondType>::IsDynamic) ||
      (TypeTrait<FirstType>::IsSparse && TypeTrait<SecondType>::IsSparse)))
    };
  };

  /**
   * Retrieve the type of a the trait of a given type (fallback).
   * \tparam T The type to be specialized.
   * \tparam N The index of the type to retrieve.
   */
  template<typename T, int N = 0>
  struct RetriveType;


  /**
   * Retrieve the type of a the trait of a given type.
   * \tparam BaseType The base template type.
   * \tparam FirstType The first template type.
   * \tparam OtherTypes Other template types.
   */
  template<template<typename...> typename BaseType, typename FirstType, typename... OtherTypes>
  struct RetriveType<BaseType<FirstType, OtherTypes...>>
  {
    using Full  = BaseType<FirstType, OtherTypes...>;
    using First = FirstType;
  };

  /*\
   |   ____  _          __  __
   |  / ___|| |_ _   _ / _|/ _|
   |  \___ \| __| | | | |_| |_
   |   ___) | |_| |_| |  _|  _|
   |  |____/ \__|\__,_|_| |_|
   |
  \*/

  /**
  * \brief Retrieve the Unicode character for the top-left corner of a table.
  * \return Unicode character for the top-left corner of a table.
  */
  static constexpr std::string table_top_left_corner() {return std::string("┌");}

  /**
  * \brief Retrieve the Unicode character for the top-right corner of a table.
  * \return Unicode character for the top-right corner of a table.
  */
  static constexpr std::string table_top_right_corner() {return std::string("┐");}

  /**
  * \brief Retrieve the Unicode character for the bottom-left corner of a table.
  * \return Unicode character for the bottom-left corner of a table.
  */
  static constexpr std::string table_bottom_left_corner() {return std::string("└");}

  /**
  * \brief Retrieve the Unicode character for the bottom-right corner of a table.
  * \return Unicode character for the bottom-right corner of a table.
  */
  static constexpr std::string table_bottom_right_corner() {return std::string("┘");}

  /**
  * \brief Retrieve the Unicode character for the left junction of a table.
  * \return Unicode character for the left junction of a table.
  */
  static constexpr std::string table_left_junction() {return std::string("├");}

  /**
  * \brief Retrieve the Unicode character for the right junction of a table.
  * \return Unicode character for the right junction of a table.
  */
  static constexpr std::string table_right_junction() {return std::string("┤");}

  /**
  * \brief Retrieve the Unicode character for the top junction of a table.
  * \return Unicode character for the top junction of a table.
  */
  static constexpr std::string table_top_junction() {return std::string("┬");}

  /**
  * \brief Retrieve the Unicode character for the bottom junction of a table.
  * \return Unicode character for the bottom junction of a table.
  */
  static constexpr std::string table_bottom_junction() {return std::string("┴");}

  /**
  * \brief Retrieve the Unicode character for the center cross of a table.
  * \return Unicode character for the center cross of a table.
  */
  static constexpr std::string table_center_cross() {return std::string("┼");}

  /**
  * \brief Retrieve the Unicode character for the horizontal line of a table.
  * \return Unicode character for the horizontal line of a table.
  */
  static constexpr std::string table_horizontal_line() {return std::string("─");}

  /**
  * \brief Retrieve the Unicode character for a number of horizontal lines of a table.
  * \return Unicode character for the horizontal lines of a table.
  * \tparam N Number of horizontal lines.
  */
  template <Integer N>
  static constexpr std::string table_horizontal_line() {
    std::string line;
    for (Integer i{0}; i < N; ++i) {line += table_horizontal_line();}
    return line;
  }

  /**
  * \brief Retrieve the Unicode character for the vertical line of a table.
  * \return Unicode character for the vertical line of a table.
  */
  static constexpr std::string table_vertical_line() {return std::string("│");}

} // namespace Optimist

#endif // INCLUDE_OPTIMIST_HH
