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

#ifndef OPTIMIST_VECTOR_FUNCTION_BOOTH_HXX
#define OPTIMIST_VECTOR_FUNCTION_BOOTH_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the Booth function.
    *
    * Class container for the Booth function, which is defined as:
    * \f[
    * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 + 2x_2 - 7 \\ 2x_1 + x_2 - 5 \end{bmatrix} \text{.}
    * \f]
    * The function has one solution at \f$\mathbf{x} = [1, 3]^\top\f$, with \f$f(\mathbf{x}) = 0\f$.
    * The initial guesses are generated on the square \f$x_i \in [-10, 10]\f$, for all \f$x_i = 1, 2\f$.
    */
    class Booth : public VectorFunction<2, 2, Booth>
    {
    public:
      using Vector = typename VectorFunction<2, 2, Booth>::InputVector;
      using Matrix = typename VectorFunction<2, 2, Booth>::Matrix;
      using Tensor = typename VectorFunction<2, 2, Booth>::Tensor;

      /**
      * Class constructor for the Booth function.
      */
      Booth()
      {
        this->m_solutions.emplace_back(1.0, 3.0);
        for (Real x{-10.0}; x < 10.0 + EPSILON; x += 5.0) {
          for (Real y{-10.0}; y < 10.0 + EPSILON; y += 5.0) {
            this->m_guesses.emplace_back(x, y);
          }
        }
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Booth";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const Vector & x, Vector & out) const
      {
        out << x(0) + 2.0*x(1) - 7.0, 2.0*x(0) + x(1) - 5.0;
      }
      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const Vector & /*x*/, Matrix & out) const
      {
        out << 1.0, 2.0, 2.0, 1.0;
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const Vector & /*x*/, Tensor & out) const
      {
        out.resize(this->output_dimension());
        std::for_each(out.begin(), out.end(), [] (Matrix& m) {m.setZero();});
      }

    }; // class Booth

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_VECTOR_FUNCTION_BOOTH_HXX
