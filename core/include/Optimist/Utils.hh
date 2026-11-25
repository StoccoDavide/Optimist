/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_UTILS_HH
#define OPTIMIST_UTILS_HH

#include "Optimist.hh"

namespace Optimist {

  /**
   * \brief Helpers to move data between the input types accepted by the
   * solvers (fixed-size, dynamic-size and sparse Eigen vectors) and the dense
   * dynamic-size storage used internally by some algorithms.
   */
  namespace Utils {

    /**
     * Dense dynamic-size column vector with the same scalar type as \c Vector.
     */
    template <typename Vector>
    using DenseVector =
        Eigen::Matrix<typename TypeTrait<Vector>::Scalar, Eigen::Dynamic, 1>;

    /**
     * Dense dynamic-size matrix with the same scalar type as \c Vector.
     */
    template <typename Vector>
    using DenseMatrix = Eigen::Matrix<typename TypeTrait<Vector>::Scalar,
                                      Eigen::Dynamic,
                                      Eigen::Dynamic>;

    /**
     * Convert any supported vector to a dense dynamic-size vector.
     * \param[in] x Input vector.
     * \return The dense copy of the vector.
     */
    template <typename Vector>
    DenseVector<Vector> to_dense(const Vector &x) {
      return DenseVector<Vector>(x);
    }

    /**
     * Convert a dense vector to any supported vector type.
     * \tparam Vector Output vector type.
     * \param[in] x Dense input vector.
     * \return The converted vector.
     */
    template <typename Vector, typename Derived>
    Vector from_dense(const Eigen::MatrixBase<Derived> &x) {
      if constexpr (TypeTrait<Vector>::IsSparse) {
        Vector out(x.size());
        out.reserve(x.size());
        for (Eigen::Index i{0}; i < x.size(); ++i) {
          out.coeffRef(i) = x(i);
        }
        return out;
      } else {
        return Vector(x);
      }
    }

    /**
     * Convert a dense matrix to a dense or sparse Eigen matrix type.
     * \tparam Matrix Output matrix type.
     * \param[in] A Dense input matrix.
     * \return The converted matrix.
     */
    template <typename Matrix, typename Derived>
    Matrix from_dense_matrix(const Eigen::MatrixBase<Derived> &A) {
      if constexpr (TypeTrait<Matrix>::IsSparse) {
        return Matrix(A.sparseView());
      } else {
        return Matrix(A);
      }
    }

    /**
     * Read a coefficient of any supported vector.
     * \param[in] x Input vector.
     * \param[in] i Index.
     * \return The \f$ i \f$-th coefficient.
     */
    template <typename Vector>
    typename TypeTrait<Vector>::Scalar coeff(const Vector &x, const Integer i) {
      if constexpr (TypeTrait<Vector>::IsSparse) {
        return x.coeff(i);
      } else {
        return x(i);
      }
    }

  }  // namespace Utils

}  // namespace Optimist

#endif  // OPTIMIST_UTILS_HH
