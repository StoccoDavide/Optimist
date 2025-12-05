/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_ROSENBROCK_HH
#define OPTIMIST_TESTSET_ROSENBROCK_HH

#include "Optimist/Function.hh"

namespace Optimist
{

  namespace TestSet
  {

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
     * The function has one solution at \f$\mathbf{x} = [1, \dots 1]^\top\f$, with \f$f(\mathbf{x}) = 0\f$.
     * The initial guess is \f$x_i = [-1.2, 1, -1.2, 1, \dots, -1.2, 1]^\top\f$.
     * \tparam Vector Eigen vector type.
     * \tparam N Input dimension (must be even).
     */
    template <typename Vector, Integer N>
    requires (N % 2 == 0) && TypeTrait<Vector>::IsEigen &&
      (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == N)
    class Rosenbrock : public Function<Vector, Vector, Rosenbrock<Vector, N>>
    {
    public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar = typename Vector::Scalar;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::Input;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::Output;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::Matrix;
      using typename Function<Vector, Vector, Rosenbrock<Vector, N>>::Tensor;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the extended Rosenbrock function.
       * \param[in] n Input dimension (must be even).
       */
      Rosenbrock()
      {
        this->m_solutions.emplace_back(Output::Ones());
        this->m_guesses.emplace_back(Input::Ones());
        for (Integer i{0}; i < N; i += 2) {
          this->m_guesses[0](i) = -1.2;
          this->m_guesses[0](i+1) = 1.0;
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {return "Rosenbrock<" + std::to_string(N) + ">";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Input const & x, Output & out) const
      {
        for (Integer i{0}; i < N; i += 2) {
          out(i)   = 10.0*(x(i+1) - x(i)*x(i));
          out(i+1) = 1.0 - x(i);
        }
        return out.allFinite();
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Input const & x, Matrix & out) const
      {
        out.setZero();
        for (Integer i{0}; i < N; i += 2) {
          out(i, i)   = -20.0*x(i);
          out(i, i+1) = 10.0;
          out(i+1, i) = -1.0;
        }
        return out.allFinite();
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Input const & /*x*/, Tensor & out) const
      {
        out.resize(this->output_dimension());
        for (size_t i{0}; i < out.size(); ++i) {
          out[i].setZero();
          for (Integer j{0}; j < N; j += 2) {
            out[i](j, j) = -20.0;
          }
        }
        return true;
      }

    }; // class Rosenbrock

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


  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_ROSENBROCK_HH
