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

#ifndef OPTIMIST_TESTSET_RASTRIGIN_HH
#define OPTIMIST_TESTSET_RASTRIGIN_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Rastrigin's function.
     *
     * \f[
     *   f(\mathbf{x}) = 10N + \sum_{i=1}^{N} \left( x_i^2 - 10\cos(2\pi x_i)
     *   \right) \text{,}
     * \f]
     * a highly multimodal function with a regular lattice of local minima
     * and global minimum \f$ f(\mathbf{0}) = 0 \f$. The usual search domain is
     * \f$ [-5.12, 5.12]^N \f$. It is meant for global optimizers.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension.
     */
    template <typename Vector, Integer N>
    class Rastrigin
        : public DenseCostFunction<Vector, N, Rastrigin<Vector, N>> {
     public:
      using Base = DenseCostFunction<Vector, N, Rastrigin<Vector, N>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Rastrigin's function.
       */
      Rastrigin() {
        this->add_solution(Dense::Zero(N));
        // Off-diagonal guesses: on the diagonal the gradient points straight
        // at the global minimizer
        Dense x0(N), x1(N);
        for (Integer i{0}; i < N; ++i) {
          x0(i) = i % 2 == 0 ? 3.1 : -2.3;
          x1(i) = i % 2 == 0 ? -4.2 : 1.7;
        }
        this->add_guess(x0);
        this->add_guess(x1);
      }

      constexpr std::string name_impl() const {
        return "Rastrigin<" + std::to_string(N) + ">";
      }

      /**
       * Usual half-width of the search domain.
       */
      static constexpr Scalar domain() {
        return 5.12;
      }

      Scalar value(const Dense &x) const {
        const Scalar w{2 * static_cast<Scalar>(M_PI)};
        return 10 * N + (x.array().square() - 10 * (w * x.array()).cos()).sum();
      }

      void gradient(const Dense &x, Dense &g) const {
        const Scalar w{2 * static_cast<Scalar>(M_PI)};
        g = 2 * x.array() + 10 * w * (w * x.array()).sin();
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        const Scalar w{2 * static_cast<Scalar>(M_PI)};
        H.setZero();
        H.diagonal() = 2 + 10 * w * w * (w * x.array()).cos();
      }

    };  // class Rastrigin

    template <typename Vector>
    using Rastrigin2 = Rastrigin<Vector, 2>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_RASTRIGIN_HH
