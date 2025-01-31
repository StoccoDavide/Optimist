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
* \brief The namespace for the Optimist library.
*
* The namespace contains all the classes and functions of the Optimist library.
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

  static Real const EPSILON        = std::numeric_limits<Real>::epsilon();            /**< Machine epsilon epsilon static constant value. */
  static Real const SQRT_EPSILON   = std::sqrt(EPSILON);                              /**< Square root of machine epsilon epsilon static constant value. */
  static Real const CBRT_EPSILON   = std::cbrt(EPSILON);                              /**< Cubic root of machine epsilon epsilon static constant value. */
  static Real const EPSILON_HIGH   = Real(1.0e-12);                                   /**< High precision epsilon static constant value. */
  static Real const EPSILON_MEDIUM = Real(1.0e-10);                                   /**< Medium precision epsilon static constant value. */
  static Real const EPSILON_LOW    = Real(1.0e-08);                                   /**< Low precision epsilon static constant value. */
  static Real const INFTY          = std::numeric_limits<Real>::infinity();           /**< Infinity static constant value. */
  static Real const QUIET_NAN      = std::numeric_limits<Real>::quiet_NaN();          /**< Not-a-number static constant value. */

  static Vector2 const NAN_VEC2      = Vector2::Constant(QUIET_NAN); /**< Not-a-number \f$ 2 \times 1 \f$ vector static constant object. */
  static Matrix2 const NAN_MAT2      = Matrix2::Constant(QUIET_NAN); /**< Not-a-number \f$ 2 \times 2 \f$ matrix static constant object. */
  static Vector2 const ZEROS_VEC2    = Vector2::Zero();              /**< Zeros \f$ 2 \times 1 \f$ vector static constant object. */
  static Matrix2 const ZEROS_MAT2    = Matrix2::Zero();              /**< Zeros \f$ 2 \times 2 \f$ matrix static constant object. */
  static Vector2 const ONES_VEC2     = Vector2::Ones();              /**< Ones \f$ 2 \times 1 \f$ vector static constant object. */
  static Matrix2 const ONES_MAT2     = Matrix2::Ones();              /**< Ones \f$ 2 \times 2 \f$ matrix static constant object. */
  static Matrix2 const IDENTITY_MAT2 = Matrix2::Identity();          /**< Identity \f$ 2 \times 2 \f$ matrix static constant object. */

  static Vector3 const NAN_VEC3      = Vector3::Constant(QUIET_NAN); /**< Not-a-number \f$ 3 \times 1 \f$ vector static constant object. */
  static Matrix3 const NAN_MAT3      = Matrix3::Constant(QUIET_NAN); /**< Not-a-number \f$ 3 \times 3 \f$ matrix static constant object. */
  static Vector3 const ZEROS_VEC3    = Vector3::Zero();              /**< Zeros \f$ 3 \times 1 \f$ vector static constant object. */
  static Matrix3 const ZEROS_MAT3    = Matrix3::Zero();              /**< Zeros \f$ 3 \times 3 \f$ matrix static constant object. */
  static Vector3 const ONES_VEC3     = Vector3::Ones();              /**< Ones \f$ 3 \times 1 \f$ vector static constant object. */
  static Matrix3 const ONES_MAT3     = Matrix3::Ones();              /**< Ones \f$ 3 \times 3 \f$ matrix static constant object. */
  static Matrix3 const IDENTITY_MAT3 = Matrix3::Identity();          /**< Identity \f$ 3 \times 3 \f$ matrix static constant object. */

  static Vector4 const NAN_VEC4      = Vector4::Constant(QUIET_NAN); /**< Not-a-number \f$ 4 \times 1 \f$ vector static constant object. */
  static Matrix4 const NAN_MAT4      = Matrix4::Constant(QUIET_NAN); /**< Not-a-number \f$ 4 \times 4 \f$ matrix static constant object. */
  static Vector4 const ZEROS_VEC4    = Vector4::Zero();              /**< Zeros \f$ 4 \times 1 \f$ vector static constant object. */
  static Matrix4 const ZEROS_MAT4    = Matrix4::Zero();              /**< Zeros \f$ 4 \times 4 \f$ matrix static constant object. */
  static Vector4 const ONES_VEC4     = Vector4::Ones();              /**< Ones \f$ 4 \times 1 \f$ vector static constant object. */
  static Matrix4 const ONES_MAT4     = Matrix4::Ones();              /**< Ones \f$ 4 \times 4 \f$ matrix static constant object. */
  static Matrix4 const IDENTITY_MAT4 = Matrix4::Identity();          /**< Identity \f$ 4 \times 4 \f$ matrix static constant object. */

  static Vector5 const NAN_VEC5      = Vector5::Constant(QUIET_NAN); /**< Not-a-number \f$ 5 \times 1 \f$ vector static constant object. */
  static Matrix5 const NAN_MAT5      = Matrix5::Constant(QUIET_NAN); /**< Not-a-number \f$ 5 \times 5 \f$ matrix static constant object. */
  static Vector5 const ZEROS_VEC5    = Vector5::Zero();              /**< Zeros \f$ 5 \times 1 \f$ vector static constant object. */
  static Matrix5 const ZEROS_MAT5    = Matrix5::Zero();              /**< Zeros \f$ 5 \times 5 \f$ matrix static constant object. */
  static Vector5 const ONES_VEC5     = Vector5::Ones();              /**< Ones \f$ 5 \times 1 \f$ vector static constant object. */
  static Matrix5 const ONES_MAT5     = Matrix5::Ones();              /**< Ones \f$ 5 \times 5 \f$ matrix static constant object. */
  static Matrix5 const IDENTITY_MAT5 = Matrix5::Identity();          /**< Identity \f$ 5 \times 5 \f$ matrix static constant object. */

  static Vector6 const NAN_VEC6      = Vector6::Constant(QUIET_NAN); /**< Not-a-number \f$ 6 \times 1 \f$ vector static constant object. */
  static Matrix6 const NAN_MAT6      = Matrix6::Constant(QUIET_NAN); /**< Not-a-number \f$ 6 \times 6 \f$ matrix static constant object. */
  static Vector6 const ZEROS_VEC6    = Vector6::Zero();              /**< Zeros \f$ 6 \times 1 \f$ vector static constant object. */
  static Matrix6 const ZEROS_MAT6    = Matrix6::Zero();              /**< Zeros \f$ 6 \times 6 \f$ matrix static constant object. */
  static Vector6 const ONES_VEC6     = Vector6::Ones();              /**< Ones \f$ 6 \times 1 \f$ vector static constant object. */
  static Matrix6 const ONES_MAT6     = Matrix6::Ones();              /**< Ones \f$ 6 \times 6 \f$ matrix static constant object. */
  static Matrix6 const IDENTITY_MAT6 = Matrix6::Identity();          /**< Identity \f$ 6 \times 6 \f$ matrix static constant object. */

  static Vector7 const NAN_VEC7      = Vector7::Constant(QUIET_NAN); /**< Not-a-number \f$ 7 \times 1 \f$ vector static constant object. */
  static Matrix7 const NAN_MAT7      = Matrix7::Constant(QUIET_NAN); /**< Not-a-number \f$ 7 \times 7 \f$ matrix static constant object. */
  static Vector7 const ZEROS_VEC7    = Vector7::Zero();              /**< Zeros \f$ 7 \times 1 \f$ vector static constant object. */
  static Matrix7 const ZEROS_MAT7    = Matrix7::Zero();              /**< Zeros \f$ 7 \times 7 \f$ matrix static constant object. */
  static Vector7 const ONES_VEC7     = Vector7::Ones();              /**< Ones \f$ 7 \times 1 \f$ vector static constant object. */
  static Matrix7 const ONES_MAT7     = Matrix7::Ones();              /**< Ones \f$ 7 \times 7 \f$ matrix static constant object. */
  static Matrix7 const IDENTITY_MAT7 = Matrix7::Identity();          /**< Identity \f$ 7 \times 7 \f$ matrix static constant object. */

  static Vector8 const NAN_VEC8      = Vector8::Constant(QUIET_NAN); /**< Not-a-number \f$ 8 \times 1 \f$ vector static constant object. */
  static Matrix8 const NAN_MAT8      = Matrix8::Constant(QUIET_NAN); /**< Not-a-number \f$ 8 \times 8 \f$ matrix static constant object. */
  static Vector8 const ZEROS_VEC8    = Vector8::Zero();              /**< Zeros \f$ 8 \times 1 \f$ vector static constant object. */
  static Matrix8 const ZEROS_MAT8    = Matrix8::Zero();              /**< Zeros \f$ 8 \times 8 \f$ matrix static constant object. */
  static Vector8 const ONES_VEC8     = Vector8::Ones();              /**< Ones \f$ 8 \times 1 \f$ vector static constant object. */
  static Matrix8 const ONES_MAT8     = Matrix8::Ones();              /**< Ones \f$ 8 \times 8 \f$ matrix static constant object. */
  static Matrix8 const IDENTITY_MAT8 = Matrix8::Identity();          /**< Identity \f$ 8 \times 8 \f$ matrix static constant object. */

  static Vector9 const NAN_VEC9      = Vector9::Constant(QUIET_NAN); /**< Not-a-number \f$ 9 \times 1 \f$ vector static constant object. */
  static Matrix9 const NAN_MAT9      = Matrix9::Constant(QUIET_NAN); /**< Not-a-number \f$ 9 \times 9 \f$ matrix static constant object. */
  static Vector9 const ZEROS_VEC9    = Vector9::Zero();              /**< Zeros \f$ 9 \times 1 \f$ vector static constant object. */
  static Matrix9 const ZEROS_MAT9    = Matrix9::Zero();              /**< Zeros \f$ 9 \times 9 \f$ matrix static constant object. */
  static Vector9 const ONES_VEC9     = Vector9::Ones();              /**< Ones \f$ 9 \times 1 \f$ vector static constant object. */
  static Matrix9 const ONES_MAT9     = Matrix9::Ones();              /**< Ones \f$ 9 \times 9 \f$ matrix static constant object. */
  static Matrix9 const IDENTITY_MAT9 = Matrix9::Identity();          /**< Identity \f$ 9 \times 9 \f$ matrix static constant object. */

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

// Time measurement
#include "Optimist/TicToc.hxx"

// The generic solver
#include "Optimist/Solver.hxx"

// Scalar solvers
#include "Optimist/ScalarSolver.hxx"
#include "Optimist/ScalarSolver/Newton.hxx"

// Root-finding solvers
#include "Optimist/RootFinder.hxx"
#include "Optimist/RootFinder/Newton.hxx"
#include "Optimist/RootFinder/Broyden.hxx"
#include "Optimist/RootFinder/Greenstadt.hxx"

// Optimization solvers
#include "Optimist/Optimizer.hxx"

// Mapping solvers to functions
#include "Optimist/Factory.hxx"

#endif // INCLUDE_OPTIMIST_HH
