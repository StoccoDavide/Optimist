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

#ifndef OPTIMIST_TESTSET_ROSENBROCK_HH
#define OPTIMIST_TESTSET_ROSENBROCK_HH

#include "Optimist/Function.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the extended Rosenbrock function.
     *
     * Class container for the extended Rosenbrock function, which defined as:
     * \f[
     * \mathbf{f}(\mathbf{x}) = \begin{bmatrix}
     *   10(x_2 - x_1^2) \\ 1 - x_1 \\
     *   10(x_4 - x_3^2) \\ 1 - x_3 \\
     *   \vdots \\
     *   10(x_N - x_{N-1}^2) \\ 1 - x_{N-1}
     * \end{bmatrix} \text{.}
     * \f]
     * The function has one solution at \f$\mathbf{x} = [1, \dots 1]^\top\f$,
     * with
     * \f$f(\mathbf{x}) = 0\f$. The initial guess is \f$x_i = [-1.2, 1, -1.2, 1,
     * \dots, -1.2, 1]^\top\f$.
     * \tparam Vector Eigen vector type.
     * \tparam N Input dimension (must be even).
     */
    template <typename Vector, Integer N>
      requires(N % 2 == 0) && TypeTrait<Vector>::IsEigen &&
              (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == N)
    class Rosenbrock : public Function<Vector, Vector, Rosenbrock<Vector, N>> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::
          FirstDerivative;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::
          SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the extended Rosenbrock function.
       */
      Rosenbrock() {
        this->m_solutions.resize(1);
        this->m_guesses.resize(1);
        if constexpr (VectorTrait::IsFixed) {
          this->m_solutions[0].setConstant(1.0);
          for (Integer i{0}; i < N; i += 2) {
            this->m_guesses[0](i)     = -1.2;
            this->m_guesses[0](i + 1) = 1.0;
          }
        } else if constexpr (VectorTrait::IsDynamic) {
          this->m_solutions[0].resize(N);
          this->m_solutions[0].setConstant(1.0);
          this->m_guesses[0].resize(N);
          for (Integer i{0}; i < N; i += 2) {
            this->m_guesses[0](i)     = -1.2;
            this->m_guesses[0](i + 1) = 1.0;
          }
        } else if constexpr (VectorTrait::IsSparse) {
          this->m_solutions[0].resize(N);
          this->m_solutions[0].reserve(N);
          for (Integer i{0}; i < N; ++i) {
            this->m_solutions[0].coeffRef(i) = 1.0;
          }
          this->m_guesses[0].resize(N);
          this->m_guesses[0].reserve(N);
          for (Integer i{0}; i < N; i += 2) {
            this->m_guesses[0].coeffRef(i)     = -1.2;
            this->m_guesses[0].coeffRef(i + 1) = 1.0;
          }
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Rosenbrock<" + std::to_string(N) + ">";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(const Vector &x, Vector &out) const {
#define CMD                                              \
  "Optimist::TestSet::Rosenbrock<" + std::to_string(N) + \
      ">::evaluate_impl(...): "

        if constexpr (VectorTrait::IsFixed) {
          for (Integer i{0}; i < N; i += 2) {
            out(i)     = 10.0 * (x(i + 1) - x(i) * x(i));
            out(i + 1) = 1.0 - x(i);
          }
          return out.allFinite();
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(N);
          for (Integer i{0}; i < N; i += 2) {
            out(i)     = 10.0 * (x(i + 1) - x(i) * x(i));
            out(i + 1) = 1.0 - x(i);
          }
          return out.allFinite();
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(N);
          out.reserve(N);
          for (Integer i{0}; i < N; i += 2) {
            out.coeffRef(i) = 10.0 * (x.coeff(i + 1) - x.coeff(i) * x.coeff(i));
            out.coeffRef(i + 1) = 1.0 - x.coeff(i);
          }
          for (typename Vector::InnerIterator it(out); it; ++it) {
            if (!std::isfinite(it.value())) {
              return false;
            }
          }
          return true;
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }

#undef CMD
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(const Vector &x, FirstDerivative &out) const {
#define CMD                                              \
  "Optimist::TestSet::Rosenbrock<" + std::to_string(N) + \
      ">::first_derivative_impl(...): "

        if constexpr (VectorTrait::IsFixed) {
          out.setZero();
          for (Integer i{0}; i < N; i += 2) {
            out(i, i)     = -20.0 * x(i);
            out(i, i + 1) = 10.0;
            out(i + 1, i) = -1.0;
          }
          return out.allFinite();
        } else if constexpr (VectorTrait::IsDynamic) {
          out.resize(N, N);
          out.setZero();
          for (Integer i{0}; i < N; i += 2) {
            out(i, i)     = -20.0 * x(i);
            out(i, i + 1) = 10.0;
            out(i + 1, i) = -1.0;
          }
          return out.allFinite();
        } else if constexpr (VectorTrait::IsSparse) {
          out.resize(N, N);
          out.reserve(N * N);
          for (Integer i{0}; i < N; i += 2) {
            out.coeffRef(i, i)     = -20.0 * x.coeff(i);
            out.coeffRef(i, i + 1) = 10.0;
            out.coeffRef(i + 1, i) = -1.0;
          }
          for (Integer k{0}; k < out.outerSize(); ++k) {
            for (typename FirstDerivative::InnerIterator it(out, k); it; ++it) {
              if (!std::isfinite(it.value())) {
                return false;
              }
            }
          }
          return true;
        } else {
          OPTIMIST_ERROR(CMD "input type not supported.");
          return false;
        }

#undef CMD
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(const Vector & /*x*/,
                                  SecondDerivative &out) const {
#define CMD                                              \
  "Optimist::TestSet::Rosenbrock<" + std::to_string(N) + \
      ">::second_derivative_impl(...): "

        out.resize(N);
        std::for_each(out.begin(), out.end(), [](auto &m) {
          if constexpr (VectorTrait::IsFixed) {
            m.setZero();
            for (Integer j{0}; j < N; j += 2) {
              m(j, j) = -20.0;
            }
          } else if constexpr (VectorTrait::IsDynamic) {
            m.resize(N, N);
            m.setZero();
            for (Integer j{0}; j < N; j += 2) {
              m(j, j) = -20.0;
            }
          } else if constexpr (VectorTrait::IsSparse) {
            m.resize(N, N);
            m.setZero();
            for (Integer j{0}; j < N; j += 2) {
              m.coeffRef(j, j) = -20.0;
            }
          } else {
            OPTIMIST_ERROR(CMD "input type not supported.");
          }
        });
        return true;

#undef CMD
      }

    };  // class Rosenbrock

    /**
     * \brief Class container for the 2D Rosenbrock function.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    using Rosenbrock2 = Rosenbrock<Vector, 2>;

    /**
     * \brief Class container for the 4D Rosenbrock function.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    using Rosenbrock4 = Rosenbrock<Vector, 4>;

    /**
     * \brief Class container for the 6D Rosenbrock function.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    using Rosenbrock6 = Rosenbrock<Vector, 6>;

    /**
     * \brief Class container for the 8D Rosenbrock function.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    using Rosenbrock8 = Rosenbrock<Vector, 8>;

    /**
     * \brief Class container for the 10D Rosenbrock function.
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    using Rosenbrock10 = Rosenbrock<Vector, 10>;

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ROSENBROCK_HH
