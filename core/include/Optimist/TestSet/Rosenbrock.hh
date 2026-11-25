

#pragma once

#ifndef OPTIMIST_TESTSET_ROSENBROCK_HH
#define OPTIMIST_TESTSET_ROSENBROCK_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the extended Rosenbrock function.
     *
     * Square system of \f$ N \f$ (even) equations,
     * \f[
     *   \mathbf{f}(\mathbf{x}) = \begin{bmatrix}
     *     10(x_2 - x_1^2) \\ 1 - x_1 \\
     *     \vdots \\
     *     10(x_N - x_{N-1}^2) \\ 1 - x_{N-1}
     *   \end{bmatrix} \text{,}
     * \f]
     * with solution \f$ \mathbf{x} = (1, \dots, 1) \f$, from the guess \f$
     * (-1.2, 1, \dots, -1.2, 1) \f$. See RosenbrockCost for the scalar cost
     * function.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension (even).
     */
    template <typename Vector, Integer N>
      requires(N % 2 == 0)
    class Rosenbrock
        : public DenseSystemFunction<Vector, N, Rosenbrock<Vector, N>> {
     public:
      using Base = DenseSystemFunction<Vector, N, Rosenbrock<Vector, N>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for the extended Rosenbrock function.
       */
      Rosenbrock() {
        Dense x0(N);
        for (Integer i{0}; i < N; i += 2) {
          x0(i)     = -1.2;
          x0(i + 1) = 1.0;
        }
        this->add_solution(Dense::Ones(N));
        this->add_guess(x0);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Rosenbrock<" + std::to_string(N) + ">";
      }

      void residual(const Dense &x, Dense &F) const {
        for (Integer i{0}; i < N; i += 2) {
          F(i)     = 10 * (x(i + 1) - x(i) * x(i));
          F(i + 1) = 1 - x(i);
        }
      }

      void jacobian(const Dense &x, DenseMatrix &J) const {
        J.setZero();
        for (Integer i{0}; i < N; i += 2) {
          J(i, i)     = -20 * x(i);
          J(i, i + 1) = 10;
          J(i + 1, i) = -1;
        }
      }

      void hessian([[maybe_unused]] const Dense &x,
                   const Integer i,
                   DenseMatrix &H) const {
        H.setZero();
        if (i % 2 == 0) {
          H(i, i) = -20;
        }
      }

    };  // class Rosenbrock

    template <typename Vector>
    using Rosenbrock2 = Rosenbrock<Vector, 2>;

    template <typename Vector>
    using Rosenbrock4 = Rosenbrock<Vector, 4>;

    template <typename Vector>
    using Rosenbrock6 = Rosenbrock<Vector, 6>;

    template <typename Vector>
    using Rosenbrock8 = Rosenbrock<Vector, 8>;

    template <typename Vector>
    using Rosenbrock10 = Rosenbrock<Vector, 10>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ROSENBROCK_HH
