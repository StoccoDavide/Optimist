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

#ifndef OPTIMIST_TESTSET_ROSENBROCK_COST_HH
#define OPTIMIST_TESTSET_ROSENBROCK_COST_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Rosenbrock's cost function in \f$ N \f$ dimensions (chained
     * form).
     *
     * \f[
     *   f(\mathbf{x}) = \sum_{i=1}^{N-1} 100(x_{i+1} - x_i^2)^2 + (1 - x_i)^2
     *   \text{,}
     * \f]
     * with global minimum \f$ f(1, \dots, 1) = 0 \f$ at the end of a long,
     * narrow, curved valley. Unlike Rosenbrock, which exposes the residuals of
     * the decoupled variant to root finders, this is the scalar cost function
     * used to benchmark optimizers.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension.
     */
    template <typename Vector, Integer N>
      requires(N >= 2)
    class RosenbrockCost
        : public DenseCostFunction<Vector, N, RosenbrockCost<Vector, N>> {
     public:
      using Base = DenseCostFunction<Vector, N, RosenbrockCost<Vector, N>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for Rosenbrock's cost function.
       */
      RosenbrockCost() {
        this->add_solution(Dense::Ones(N));
        // The classic (-1.2, 1, ..., 1) start is only used for N = 2: for
        // 4 <= N <= 7 it leads to a genuine local minimizer near x_1 = -1
        if constexpr (N == 2) {
          this->add_guess(this->make({-1.2, 1.0}));
        } else {
          this->add_guess(Dense::Constant(N, -0.5));
        }
        this->add_guess(Dense::Zero(N));
        this->add_guess(Dense::Constant(N, 2.0));
      }

      constexpr std::string name_impl() const {
        return "RosenbrockCost<" + std::to_string(N) + ">";
      }

      Scalar value(const Dense &x) const {
        Scalar f{0.0};
        for (Integer i{0}; i < N - 1; ++i) {
          const Scalar a{x(i + 1) - x(i) * x(i)}, b{1 - x(i)};
          f += 100 * a * a + b * b;
        }
        return f;
      }

      void gradient(const Dense &x, Dense &g) const {
        g.setZero();
        for (Integer i{0}; i < N - 1; ++i) {
          const Scalar a{x(i + 1) - x(i) * x(i)};
          g(i) += -400 * x(i) * a - 2 * (1 - x(i));
          g(i + 1) += 200 * a;
        }
      }

      void hessian(const Dense &x, DenseMatrix &H) const {
        H.setZero();
        for (Integer i{0}; i < N - 1; ++i) {
          H(i, i) += 1200 * x(i) * x(i) - 400 * x(i + 1) + 2;
          H(i, i + 1) += -400 * x(i);
          H(i + 1, i) += -400 * x(i);
          H(i + 1, i + 1) += 200;
        }
      }

    };  // class RosenbrockCost

    template <typename Vector>
    using RosenbrockCost2 = RosenbrockCost<Vector, 2>;

    template <typename Vector>
    using RosenbrockCost4 = RosenbrockCost<Vector, 4>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ROSENBROCK_COST_HH
