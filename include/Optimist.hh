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
#include <chrono>

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

/**
* \brief Namespace for the Optimist library.
*/
namespace Optimist
{

  /*\
   |      _    _ _
   |     / \  | (_) __ _ ___  ___  ___
   |    / _ \ | | |/ _` / __|/ _ \/ __|
   |   / ___ \| | | (_| \__ \  __/\__ \
   |  /_/   \_\_|_|\__,_|___/\___||___/
   |
  \*/

  using Real    = double; /**< Real number type. */
  using Integer = int;    /**< Integer number type. */

  using Vector2 = Eigen::Vector<Real, 2>;    /**< \f$ 2 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix2 = Eigen::Matrix<Real, 2, 2>; /**< \f$ 2 \times 2 \f$ matrix of Real number type. */
  using Vector3 = Eigen::Vector<Real, 3>;    /**< \f$ 3 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix3 = Eigen::Matrix<Real, 3, 3>; /**< \f$ 3 \times 3 \f$ matrix of Real number type. */
  using Vector4 = Eigen::Vector<Real, 4>;    /**< \f$ 4 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix4 = Eigen::Matrix<Real, 4, 4>; /**< \f$ 4 \times 4 \f$ matrix of Real number type. */
  using Vector5 = Eigen::Vector<Real, 5>;    /**< \f$ 5 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix5 = Eigen::Matrix<Real, 5, 5>; /**< \f$ 5 \times 5 \f$ matrix of Real number type. */
  using Vector6 = Eigen::Vector<Real, 6>;    /**< \f$ 6 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix6 = Eigen::Matrix<Real, 6, 6>; /**< \f$ 6 \times 6 \f$ matrix of Real number type. */
  using Vector7 = Eigen::Vector<Real, 7>;    /**< \f$ 7 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix7 = Eigen::Matrix<Real, 7, 7>; /**< \f$ 7 \times 7 \f$ matrix of Real number type. */
  using Vector8 = Eigen::Vector<Real, 8>;    /**< \f$ 8 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix8 = Eigen::Matrix<Real, 8, 8>; /**< \f$ 8 \times 8 \f$ matrix of Real number type. */
  using Vector9 = Eigen::Vector<Real, 9>;    /**< \f$ 9 \times 1 \f$ vector of Real number type (column vector). */
  using Matrix9 = Eigen::Matrix<Real, 9, 9>; /**< \f$ 9 \times 9 \f$ matrix of Real number type. */

  /**
  * \f$ N \times 1 \f$ vector of Real number type (column vector).
  * \tparam N The dimension of the vector.
  */
  template <Integer N>
  using Vector = Eigen::Vector<Real, N>;
  /**
  * \f$ N \times M \f$ matrix of Real number type.
  * \tparam N The row dimension of the matrix.
  * \tparam M The columns dimension of the matrix.
  */
  template <Integer N, Integer M>
  using Matrix = Eigen::Matrix<Real, N, M>;

  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;                 /**< \f$ N \times 1 \f$ vector of Real number type (column vector). */
  using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; /**< \f$ N \times N \f$ matrix of Real number type. */

  /*\
   |    ____                _              _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
   |
  \*/

  static const Real EPSILON        = std::numeric_limits<Real>::epsilon();            /**< Machine epsilon epsilon static constant value. */
  static const Real SQRT_EPSILON   = std::sqrt(EPSILON);                              /**< Square root of machine epsilon epsilon static constant value. */
  static const Real CBRT_EPSILON   = std::cbrt(EPSILON);                              /**< Cubic root of machine epsilon epsilon static constant value. */
  static const Real EPSILON_HIGH   = Real(1.0e-12);                                   /**< High precision epsilon static constant value. */
  static const Real EPSILON_MEDIUM = Real(1.0e-10);                                   /**< Medium precision epsilon static constant value. */
  static const Real EPSILON_LOW    = Real(1.0e-08);                                   /**< Low precision epsilon static constant value. */
  static const Real INFTY          = std::numeric_limits<Real>::infinity();           /**< Infinity static constant value. */
  static const Real QUIET_NAN      = std::numeric_limits<Real>::quiet_NaN();          /**< Not-a-number static constant value. */

  static const Vector2 NAN_VEC2      = Vector2::Constant(QUIET_NAN); /**< Not-a-number \f$ 2 \times 1 \f$ vector static constant object. */
  static const Matrix2 NAN_MAT2      = Matrix2::Constant(QUIET_NAN); /**< Not-a-number \f$ 2 \times 2 \f$ matrix static constant object. */
  static const Vector2 ZEROS_VEC2    = Vector2::Zero();              /**< Zeros \f$ 2 \times 1 \f$ vector static constant object. */
  static const Matrix2 ZEROS_MAT2    = Matrix2::Zero();              /**< Zeros \f$ 2 \times 2 \f$ matrix static constant object. */
  static const Vector2 ONES_VEC2     = Vector2::Ones();              /**< Ones \f$ 2 \times 1 \f$ vector static constant object. */
  static const Matrix2 ONES_MAT2     = Matrix2::Ones();              /**< Ones \f$ 2 \times 2 \f$ matrix static constant object. */
  static const Matrix2 IDENTITY_MAT2 = Matrix2::Identity();          /**< Identity \f$ 2 \times 2 \f$ matrix static constant object. */

  static const Vector3 NAN_VEC3      = Vector3::Constant(QUIET_NAN); /**< Not-a-number \f$ 3 \times 1 \f$ vector static constant object. */
  static const Matrix3 NAN_MAT3      = Matrix3::Constant(QUIET_NAN); /**< Not-a-number \f$ 3 \times 3 \f$ matrix static constant object. */
  static const Vector3 ZEROS_VEC3    = Vector3::Zero();              /**< Zeros \f$ 3 \times 1 \f$ vector static constant object. */
  static const Matrix3 ZEROS_MAT3    = Matrix3::Zero();              /**< Zeros \f$ 3 \times 3 \f$ matrix static constant object. */
  static const Vector3 ONES_VEC3     = Vector3::Ones();              /**< Ones \f$ 3 \times 1 \f$ vector static constant object. */
  static const Matrix3 ONES_MAT3     = Matrix3::Ones();              /**< Ones \f$ 3 \times 3 \f$ matrix static constant object. */
  static const Matrix3 IDENTITY_MAT3 = Matrix3::Identity();          /**< Identity \f$ 3 \times 3 \f$ matrix static constant object. */

  static const Vector4 NAN_VEC4      = Vector4::Constant(QUIET_NAN); /**< Not-a-number \f$ 4 \times 1 \f$ vector static constant object. */
  static const Matrix4 NAN_MAT4      = Matrix4::Constant(QUIET_NAN); /**< Not-a-number \f$ 4 \times 4 \f$ matrix static constant object. */
  static const Vector4 ZEROS_VEC4    = Vector4::Zero();              /**< Zeros \f$ 4 \times 1 \f$ vector static constant object. */
  static const Matrix4 ZEROS_MAT4    = Matrix4::Zero();              /**< Zeros \f$ 4 \times 4 \f$ matrix static constant object. */
  static const Vector4 ONES_VEC4     = Vector4::Ones();              /**< Ones \f$ 4 \times 1 \f$ vector static constant object. */
  static const Matrix4 ONES_MAT4     = Matrix4::Ones();              /**< Ones \f$ 4 \times 4 \f$ matrix static constant object. */
  static const Matrix4 IDENTITY_MAT4 = Matrix4::Identity();          /**< Identity \f$ 4 \times 4 \f$ matrix static constant object. */

  static const Vector5 NAN_VEC5      = Vector5::Constant(QUIET_NAN); /**< Not-a-number \f$ 5 \times 1 \f$ vector static constant object. */
  static const Matrix5 NAN_MAT5      = Matrix5::Constant(QUIET_NAN); /**< Not-a-number \f$ 5 \times 5 \f$ matrix static constant object. */
  static const Vector5 ZEROS_VEC5    = Vector5::Zero();              /**< Zeros \f$ 5 \times 1 \f$ vector static constant object. */
  static const Matrix5 ZEROS_MAT5    = Matrix5::Zero();              /**< Zeros \f$ 5 \times 5 \f$ matrix static constant object. */
  static const Vector5 ONES_VEC5     = Vector5::Ones();              /**< Ones \f$ 5 \times 1 \f$ vector static constant object. */
  static const Matrix5 ONES_MAT5     = Matrix5::Ones();              /**< Ones \f$ 5 \times 5 \f$ matrix static constant object. */
  static const Matrix5 IDENTITY_MAT5 = Matrix5::Identity();          /**< Identity \f$ 5 \times 5 \f$ matrix static constant object. */

  static const Vector6 NAN_VEC6      = Vector6::Constant(QUIET_NAN); /**< Not-a-number \f$ 6 \times 1 \f$ vector static constant object. */
  static const Matrix6 NAN_MAT6      = Matrix6::Constant(QUIET_NAN); /**< Not-a-number \f$ 6 \times 6 \f$ matrix static constant object. */
  static const Vector6 ZEROS_VEC6    = Vector6::Zero();              /**< Zeros \f$ 6 \times 1 \f$ vector static constant object. */
  static const Matrix6 ZEROS_MAT6    = Matrix6::Zero();              /**< Zeros \f$ 6 \times 6 \f$ matrix static constant object. */
  static const Vector6 ONES_VEC6     = Vector6::Ones();              /**< Ones \f$ 6 \times 1 \f$ vector static constant object. */
  static const Matrix6 ONES_MAT6     = Matrix6::Ones();              /**< Ones \f$ 6 \times 6 \f$ matrix static constant object. */
  static const Matrix6 IDENTITY_MAT6 = Matrix6::Identity();          /**< Identity \f$ 6 \times 6 \f$ matrix static constant object. */

  static const Vector7 NAN_VEC7      = Vector7::Constant(QUIET_NAN); /**< Not-a-number \f$ 7 \times 1 \f$ vector static constant object. */
  static const Matrix7 NAN_MAT7      = Matrix7::Constant(QUIET_NAN); /**< Not-a-number \f$ 7 \times 7 \f$ matrix static constant object. */
  static const Vector7 ZEROS_VEC7    = Vector7::Zero();              /**< Zeros \f$ 7 \times 1 \f$ vector static constant object. */
  static const Matrix7 ZEROS_MAT7    = Matrix7::Zero();              /**< Zeros \f$ 7 \times 7 \f$ matrix static constant object. */
  static const Vector7 ONES_VEC7     = Vector7::Ones();              /**< Ones \f$ 7 \times 1 \f$ vector static constant object. */
  static const Matrix7 ONES_MAT7     = Matrix7::Ones();              /**< Ones \f$ 7 \times 7 \f$ matrix static constant object. */
  static const Matrix7 IDENTITY_MAT7 = Matrix7::Identity();          /**< Identity \f$ 7 \times 7 \f$ matrix static constant object. */

  static const Vector8 NAN_VEC8      = Vector8::Constant(QUIET_NAN); /**< Not-a-number \f$ 8 \times 1 \f$ vector static constant object. */
  static const Matrix8 NAN_MAT8      = Matrix8::Constant(QUIET_NAN); /**< Not-a-number \f$ 8 \times 8 \f$ matrix static constant object. */
  static const Vector8 ZEROS_VEC8    = Vector8::Zero();              /**< Zeros \f$ 8 \times 1 \f$ vector static constant object. */
  static const Matrix8 ZEROS_MAT8    = Matrix8::Zero();              /**< Zeros \f$ 8 \times 8 \f$ matrix static constant object. */
  static const Vector8 ONES_VEC8     = Vector8::Ones();              /**< Ones \f$ 8 \times 1 \f$ vector static constant object. */
  static const Matrix8 ONES_MAT8     = Matrix8::Ones();              /**< Ones \f$ 8 \times 8 \f$ matrix static constant object. */
  static const Matrix8 IDENTITY_MAT8 = Matrix8::Identity();          /**< Identity \f$ 8 \times 8 \f$ matrix static constant object. */

  static const Vector9 NAN_VEC9      = Vector9::Constant(QUIET_NAN); /**< Not-a-number \f$ 9 \times 1 \f$ vector static constant object. */
  static const Matrix9 NAN_MAT9      = Matrix9::Constant(QUIET_NAN); /**< Not-a-number \f$ 9 \times 9 \f$ matrix static constant object. */
  static const Vector9 ZEROS_VEC9    = Vector9::Zero();              /**< Zeros \f$ 9 \times 1 \f$ vector static constant object. */
  static const Matrix9 ZEROS_MAT9    = Matrix9::Zero();              /**< Zeros \f$ 9 \times 9 \f$ matrix static constant object. */
  static const Vector9 ONES_VEC9     = Vector9::Ones();              /**< Ones \f$ 9 \times 1 \f$ vector static constant object. */
  static const Matrix9 ONES_MAT9     = Matrix9::Ones();              /**< Ones \f$ 9 \times 9 \f$ matrix static constant object. */
  static const Matrix9 IDENTITY_MAT9 = Matrix9::Identity();          /**< Identity \f$ 9 \times 9 \f$ matrix static constant object. */

  /**
  * \brief Retrieve the Unicode character for the top-left corner of a table.
  * \return Unicode character for the top-left corner of a table.
  */
  static std::string table_top_left_corner() {return "┌";}

  /**
  * \brief Retrieve the Unicode character for the top-right corner of a table.
  * \return Unicode character for the top-right corner of a table.
  */
  static std::string table_top_right_corner() {return "┐";}

  /**
  * \brief Retrieve the Unicode character for the bottom-left corner of a table.
  * \return Unicode character for the bottom-left corner of a table.
  */
  static std::string table_bottom_left_corner() {return "└";}

  /**
  * \brief Retrieve the Unicode character for the bottom-right corner of a table.
  * \return Unicode character for the bottom-right corner of a table.
  */
  static std::string table_bottom_right_corner() {return "┘";}

  /**
  * \brief Retrieve the Unicode character for the left junction of a table.
  * \return Unicode character for the left junction of a table.
  */
  static std::string table_left_junction() {return "├";}

  /**
  * \brief Retrieve the Unicode character for the right junction of a table.
  * \return Unicode character for the right junction of a table.
  */
  static std::string table_right_junction() {return "┤";}

  /**
  * \brief Retrieve the Unicode character for the top junction of a table.
  * \return Unicode character for the top junction of a table.
  */
  static std::string table_top_junction() {return "┬";}

  /**
  * \brief Retrieve the Unicode character for the bottom junction of a table.
  * \return Unicode character for the bottom junction of a table.
  */
  static std::string table_bottom_junction() {return "┴";}

  /**
  * \brief Retrieve the Unicode character for the center cross of a table.
  * \return Unicode character for the center cross of a table.
  */
  static std::string table_center_cross() {return "┼";}

  /**
  * \brief Retrieve the Unicode character for the horizontal line of a table.
  * \return Unicode character for the horizontal line of a table.
  */
  static std::string table_horizontal_line() {return "─";}

  /**
  * \brief Retrieve the Unicode character for a number of horizontal lines of a table.
  * \return Unicode character for the horizontal lines of a table.
  * \tparam N Number of horizontal lines.
  */
  template <Integer N>
  static std::string table_horizontal_line() {
    std::string line;
    for (Integer i = 0; i < N; ++i) {line += table_horizontal_line();}
    return line;
  }

  /**
  * \brief Retrieve the Unicode character for the vertical line of a table.
  * \return Unicode character for the vertical line of a table.
  */
  static std::string table_vertical_line() {return "│";}

  /**
  * Print Optimist library information on a string.
  * \return A string with the Optimist library information.
  */
  std::string Info() {
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
  void Info(std::ostream &os) {os << Info();}

} // namespace Optimist

// Finite difference
//#include "Optimist/FiniteDifference.hxx"

// The functions
#include "Optimist/Function.hxx"
#include "Optimist/Function/ScalarFunction.hxx"
#include "Optimist/Function/VectorFunction.hxx"
#include "Optimist/Function/CostFunction.hxx"
#include "Optimist/Function/CutFunction.hxx"

// The solvers
#include "Optimist/Solver.hxx"

// Scalar root-finding solvers
#include "Optimist/ScalarRootFinder.hxx"
#include "Optimist/ScalarRootFinder/Newton.hxx"

// Root-finding solvers
#include "Optimist/RootFinder.hxx"
#include "Optimist/RootFinder/Newton.hxx"
#include "Optimist/RootFinder/Broyden.hxx"
#include "Optimist/RootFinder/Greenstadt.hxx"

// Optimization solvers
#include "Optimist/Optimizer.hxx"

// Scalar optimization solvers
#include "Optimist/ScalarOptimizer.hxx"

#endif // INCLUDE_OPTIMIST_HH
