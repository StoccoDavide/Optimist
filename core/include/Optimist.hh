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

#ifndef INCLUDE_OPTIMIST_HH
#define INCLUDE_OPTIMIST_HH

// C++ standard libraries
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

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
  if (!(COND)) {                   \
    OPTIMIST_ERROR(MSG);           \
  }
#endif

// Warning for Optimist
#ifndef OPTIMIST_WARNING
#define OPTIMIST_WARNING(MSG) \
  { std::cout << MSG << std::endl; }
#endif

// Warning assert for Optimist
#ifndef OPTIMIST_ASSERT_WARNING
#define OPTIMIST_ASSERT_WARNING(COND, MSG) \
  if (!(COND)) {                           \
    OPTIMIST_WARNING(MSG);                 \
  }
#endif

// Define the basic constants for Optimist
#ifndef OPTIMIST_BASIC_CONSTANTS
#define OPTIMIST_BASIC_CONSTANTS(Scalar)                                       \
  static constexpr Scalar EPSILON{                                             \
    std::numeric_limits<Scalar>::epsilon()};   /**< Machine epsilon epsilon    \
                                                  static constant value. */    \
  static constexpr Scalar INFTY{                                               \
    std::numeric_limits<Scalar>::infinity()};  /**< Infinity static constant   \
                                                  value. */                    \
  static constexpr Scalar QUIET_NAN{                                           \
    std::numeric_limits<Scalar>::quiet_NaN()}; /**< Not-a-number static        \
                                                  constant value. */           \
  inline static const Scalar SQRT_EPSILON{                                     \
    std::sqrt(EPSILON)}; /**< Square root of epsilon static constant value. */ \
  inline static const Scalar CBRT_EPSILON{                                     \
    std::cbrt(EPSILON)}; /**< Cube root of epsilon static constant value. */
#endif

#ifndef OPTIMIST_DEFAULT_INTEGER_TYPE
#define OPTIMIST_DEFAULT_INTEGER_TYPE int
#endif

/**
 * \brief Namespace for the Optimist library.
 */
namespace Optimist {

  /**
   * \brief The Integer type as used for the API.
   *
   * The Integer type, \c \#define the preprocessor symbol \c
   * OPTIMIST_DEFAULT_INTEGER_TYPE. The default value is \c int.
   */
  using Integer = OPTIMIST_DEFAULT_INTEGER_TYPE;

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
  struct TypeTrait<
      ScalarType,
      std::enable_if_t<std::is_floating_point<ScalarType>::value>> {
    using Scalar = ScalarType;
    using Type   = ScalarType;
    static constexpr Integer Dimension{1};
    static constexpr bool IsScalar{true};
    static constexpr bool IsEigen{false};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{false};
  };

  /**
   * Traits class for fixed-size dense Eigen matrices.
   * \tparam Scalar The scalar type.
   * \tparam Rows The number of rows.
   * \tparam Cols The number of columns.
   */
  template <typename ScalarType, Integer N, Integer M>
  struct TypeTrait<Eigen::Matrix<ScalarType, N, M>,
                   std::enable_if_t<(N > 0 && M > 0)>> {
    using Scalar = ScalarType;
    using Type   = Eigen::Matrix<ScalarType, N, M>;
    static constexpr Integer Dimension{N * M};
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
  struct TypeTrait<Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>> {
    using Scalar = ScalarType;
    using Type   = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
    static constexpr Integer Dimension{Eigen::Dynamic};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{true};
    static constexpr bool IsSparse{false};
  };

  /**
   * Traits class for dynamic-size dense Eigen matrices.
   * \tparam Scalar The scalar type.
   */
  template <typename ScalarType>
  struct TypeTrait<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> {
    using Scalar = ScalarType;
    using Type   = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
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
  struct TypeTrait<Eigen::SparseVector<ScalarType>> {
    using Scalar = ScalarType;
    using Type   = Eigen::SparseVector<ScalarType>;
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
  struct TypeTrait<Eigen::SparseMatrix<ScalarType, Options, StorageIndex>> {
    using Scalar = ScalarType;
    using Type   = Eigen::SparseMatrix<ScalarType, Options, StorageIndex>;
    static constexpr Integer Dimension{Eigen::Dynamic};
    static constexpr bool IsScalar{false};
    static constexpr bool IsEigen{true};
    static constexpr bool IsFixed{false};
    static constexpr bool IsDynamic{false};
    static constexpr bool IsSparse{true};
  };

  /**
   * Retrieve the type of a the trait of a given type (fallback).
   * \tparam T The type to be specialized.
   */

  /**
   * Retrieve the type of a the trait of a given type.
   * \tparam BaseType The base template type.
   * \tparam FirstType The first template type.
   * \tparam OtherTypes Other template types.
   */
  template <typename T>
  struct RetrieveType;  // primary template

  template <template <typename, auto...> class BaseType,
            typename FirstType,
            auto... Rest>
  struct RetrieveType<BaseType<FirstType, Rest...>> {
    using Full  = BaseType<FirstType, Rest...>;
    using First = FirstType;
  };

  /**
   * \brief Retrieve the Unicode character for the top-left corner of a table.
   * \return Unicode character for the top-left corner of a table.
   */
  static constexpr std::string_view table_top_left_corner() {
    return "┌";
  }

  /**
   * \brief Retrieve the Unicode character for the top-right corner of a table.
   * \return Unicode character for the top-right corner of a table.
   */
  static constexpr std::string_view table_top_right_corner() {
    return "┐";
  }

  /**
   * \brief Retrieve the Unicode character for the bottom-left corner of a
   * table.
   * \return Unicode character for the bottom-left corner of a table.
   */
  static constexpr std::string_view table_bottom_left_corner() {
    return "└";
  }

  /**
   * \brief Retrieve the Unicode character for the bottom-right corner of a
   * table.
   * \return Unicode character for the bottom-right corner of a table.
   */
  static constexpr std::string_view table_bottom_right_corner() {
    return "┘";
  }

  /**
   * \brief Retrieve the Unicode character for the left junction of a table.
   * \return Unicode character for the left junction of a table.
   */
  static constexpr std::string_view table_left_junction() {
    return "├";
  }

  /**
   * \brief Retrieve the Unicode character for the right junction of a table.
   * \return Unicode character for the right junction of a table.
   */
  static constexpr std::string_view table_right_junction() {
    return "┤";
  }

  /**
   * \brief Retrieve the Unicode character for the top junction of a table.
   * \return Unicode character for the top junction of a table.
   */
  static constexpr std::string_view table_top_junction() {
    return "┬";
  }

  /**
   * \brief Retrieve the Unicode character for the bottom junction of a table.
   * \return Unicode character for the bottom junction of a table.
   */
  static constexpr std::string_view table_bottom_junction() {
    return "┴";
  }

  /**
   * \brief Retrieve the Unicode character for the center cross of a table.
   * \return Unicode character for the center cross of a table.
   */
  static constexpr std::string_view table_center_cross() {
    return "┼";
  }

  /**
   * \brief Retrieve the Unicode character for the horizontal line of a table.
   * \return Unicode character for the horizontal line of a table.
   */
  static constexpr std::string_view table_horizontal_line() {
    return "─";
  }

  /**
   * \brief Retrieve the Unicode character for a number of horizontal lines of a
   * table.
   * \return Unicode character for the horizontal lines of a table.
   * \tparam N Number of horizontal lines.
   */
  template <Integer N>
  static std::string table_horizontal_line() {
    std::string line;
    for (Integer i{0}; i < N; ++i) {
      line += table_horizontal_line();
    }
    return line;
  }

  /**
   * \brief Retrieve the Unicode character for the vertical line of a table.
   * \return Unicode character for the vertical line of a table.
   */
  static constexpr std::string_view table_vertical_line() {
    return "│";
  }

}  // namespace Optimist

#endif  // INCLUDE_OPTIMIST_HH
