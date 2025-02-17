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

#ifndef OPTIMIST_VECTOR_FUNCTION_ROSENBROCK_HXX
#define OPTIMIST_VECTOR_FUNCTION_ROSENBROCK_HXX

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
    */
    template <Integer N>
    class Rosenbrock : public VectorFunction<N, N, Rosenbrock<N>>
    {
      static_assert(N > 0 && N % 2 == 0, "please use an even number of dimensions");

    public:
      using Vector = typename VectorFunction<N, N, Rosenbrock<N>>::InputVector;
      using Matrix = typename VectorFunction<N, N, Rosenbrock<N>>::Matrix;
      using Tensor = typename VectorFunction<N, N, Rosenbrock<N>>::Tensor;

      /**
      * Class constructor for the extended Rosenbrock function.
      */
      Rosenbrock()
      {
        this->m_solutions.emplace_back(Vector::Ones());
        this->m_guesses.emplace_back(Vector::Ones());
        for (Integer i{0}; i < N; i += 2) {
          this->m_guesses[0](i) = -1.2;
          this->m_guesses[0](i+1) = 1.0;
        }
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Rosenbrock<" + std::to_string(N) + ">";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const Vector & x, Vector & out) const
      {
        for (Integer i{0}; i < N; i += 2) {
          out(i)   = 10.0*(x(i+1) - x(i)*x(i));
          out(i+1) = 1.0 - x(i);
        }
      }
      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const Vector & x, Matrix & out) const
      {
        out.setZero();
        for (Integer i{0}; i < N; i += 2) {
          out(i, i)   = -20.0*x(i);
          out(i, i+1) = 10.0;
          out(i+1, i) = -1.0;
        }
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const Vector & /*x*/, Tensor & out) const
      {
        out.resize(this->output_dimension());
        for (Integer i{0}; i < static_cast<Integer>(out.size()); ++i) {
          out[i].setZero();
          for (Integer j{0}; j < N; j += 2) {
            out[i](j, j) = -20.0;
          }
        }
      }

    }; // class Rosenbrock<N>

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_VECTOR_FUNCTION_ROSENBROCK_HXX
