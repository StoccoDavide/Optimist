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

#ifndef OPTIMIST_TESTSET_TRID_HH
#define OPTIMIST_TESTSET_TRID_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Trid function.
     *
     * \f[
     *   f(\mathbf{x}) = \sum_{i=1}^{N} (x_i - 1)^2 -
     *   \sum_{i=2}^{N} x_i x_{i-1} \text{,}
     * \f]
     * a strictly convex quadratic whose minimum \f$ f^\star = -N(N+4)(N-1)/6
     * \f$, attained at \f$ x_i = i(N+1-i) \f$, is \em negative. It checks that
     * the cost function is minimized directly and not through a least-squares
     * reformulation.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension.
     */
    template <typename Vector, Integer N>
    class Trid : public DenseCostFunction<Vector, N, Trid<Vector, N>> {
     public:
      using Base = DenseCostFunction<Vector, N, Trid<Vector, N>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for the Trid function.
       */
      Trid() {
        Dense x(N);
        for (Integer i{0}; i < N; ++i) {
          x(i) = static_cast<Scalar>((i + 1) * (N - i));
        }
        this->add_solution(x);
        this->add_guess(Dense::Zero(N));
        this->add_guess(Dense::Constant(N, 10.0));
      }

      constexpr std::string name_impl() const {
        return "Trid<" + std::to_string(N) + ">";
      }

      /**
       * Minimum value of the function.
       */
      static constexpr Scalar minimum() {
        return -static_cast<Scalar>(N * (N + 4) * (N - 1)) / 6;
      }

      Scalar value(const Dense &x) const {
        Scalar f{(x.array() - 1).square().sum()};
        for (Integer i{1}; i < N; ++i) {
          f -= x(i) * x(i - 1);
        }
        return f;
      }

      void gradient(const Dense &x, Dense &g) const {
        for (Integer i{0}; i < N; ++i) {
          g(i) = 2 * (x(i) - 1);
          if (i > 0) {
            g(i) -= x(i - 1);
          }
          if (i < N - 1) {
            g(i) -= x(i + 1);
          }
        }
      }

      void hessian([[maybe_unused]] const Dense &x, DenseMatrix &H) const {
        H.setZero();
        for (Integer i{0}; i < N; ++i) {
          H(i, i) = 2;
          if (i > 0) {
            H(i, i - 1) = -1;
            H(i - 1, i) = -1;
          }
        }
      }

    };  // class Trid

    template <typename Vector>
    using Trid2 = Trid<Vector, 2>;

    template <typename Vector>
    using Trid6 = Trid<Vector, 6>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_TRID_HH
