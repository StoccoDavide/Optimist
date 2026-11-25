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

#ifndef OPTIMIST_TESTSET_ACKLEY_HH
#define OPTIMIST_TESTSET_ACKLEY_HH

#include <unsupported/Eigen/AutoDiff>

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Ackley's function.
     *
     * \f[
     *   f(\mathbf{x}) = -20\exp\left(-0.2\sqrt{\frac{1}{N}\sum x_i^2}\right) -
     *   \exp\left(\frac{1}{N}\sum\cos(2\pi x_i)\right) + 20 + e \text{,}
     * \f]
     * a nearly flat outer region riddled with local minima around a deep
     * funnel with global minimum \f$ f(\mathbf{0}) = 0 \f$. The usual search
     * domain is \f$ [-5, 5]^N \f$. It is meant for global optimizers. The
     * Hessian is obtained by forward-mode automatic differentiation of the
     * gradient (Eigen's \c AutoDiffScalar). The function is not
     * differentiable at the origin, where the gradient is set to zero.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension.
     */
    template <typename Vector, Integer N>
    class Ackley : public DenseCostFunction<Vector, N, Ackley<Vector, N>> {
     public:
      using Base = DenseCostFunction<Vector, N, Ackley<Vector, N>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Ackley's function.
       */
      Ackley() {
        this->add_solution(Dense::Zero(N));
        // Off-diagonal guesses: on the diagonal the gradient points straight
        // at the global minimizer
        Dense x0(N), x1(N);
        for (Integer i{0}; i < N; ++i) {
          x0(i) = i % 2 == 0 ? 2.5 : -1.3;
          x1(i) = i % 2 == 0 ? -3.5 : 3.2;
        }
        this->add_guess(x0);
        this->add_guess(x1);
      }

      constexpr std::string name_impl() const {
        return "Ackley<" + std::to_string(N) + ">";
      }

      /**
       * Usual half-width of the search domain.
       */
      static constexpr Scalar domain() {
        return 5.0;
      }

      Scalar value(const Dense &x) const {
        const Scalar n{static_cast<Scalar>(N)};
        const Scalar w{2 * static_cast<Scalar>(M_PI)};
        const Scalar r{std::sqrt(x.squaredNorm() / n)};
        const Scalar c{(w * x.array()).cos().sum() / n};
        return -20 * std::exp(-r / 5) - std::exp(c) + 20 +
               static_cast<Scalar>(M_E);
      }

      void gradient(const Dense &x, Dense &g) const {
        g = gradient_t(x);
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        // Seed one derivative direction per coordinate and differentiate the
        // gradient: H(:, j) = d(grad f)/dx_j
        using AD = Eigen::AutoDiffScalar<Dense>;
        Eigen::Matrix<AD, Eigen::Dynamic, 1> x_ad(N);
        for (Integer j{0}; j < N; ++j) {
          x_ad(j) = AD(x(j), N, j);
        }
        const Eigen::Matrix<AD, Eigen::Dynamic, 1> g_ad{gradient_t(x_ad)};
        for (Integer i{0}; i < N; ++i) {
          H.row(i) = g_ad(i).derivatives().transpose();
        }
        H = (0.5 * (H + H.transpose())).eval();
      }

     private:
      /**
       * Gradient written for a generic scalar type, so that it can be
       * differentiated automatically.
       */
      template <typename T>
      static Eigen::Matrix<T, Eigen::Dynamic, 1> gradient_t(
          const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) {
        using std::cos;
        using std::exp;
        using std::sin;
        using std::sqrt;
        const Scalar n{static_cast<Scalar>(N)};
        const Scalar w{2 * static_cast<Scalar>(M_PI)};
        T r2{x(0) * x(0)}, c{cos(w * x(0))};
        for (Integer i{1}; i < N; ++i) {
          r2 += x(i) * x(i);
          c += cos(w * x(i));
        }
        r2 /= n;
        c /= n;
        Eigen::Matrix<T, Eigen::Dynamic, 1> g(N);
        const T e_c{exp(c)};
        const bool origin{r2 == Scalar(0)};
        const T r{origin ? T(Scalar(1)) : T(sqrt(r2))};
        const T a{origin ? T(Scalar(0)) : T(4 * exp(-r / 5) / (n * r))};
        for (Integer i{0}; i < N; ++i) {
          g(i) = e_c * (w / n) * sin(w * x(i)) + a * x(i);
        }
        return g;
      }

    };  // class Ackley

    template <typename Vector>
    using Ackley2 = Ackley<Vector, 2>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ACKLEY_HH
