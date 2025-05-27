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

#ifndef OPTIMIST_TESTSET_BROWN_HH
#define OPTIMIST_TESTSET_BROWN_HH

#include "Optimist/TestSet.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the Brown badly scaled function.
    *
    * Class container for the Brown badly scaled function, which is defined as:
    * \f[
    * \mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 - a \\ x_2 - 2a \\ x_1x_2 - 2 \end{bmatrix} \text{,}
    * \f]
    * where \f$a = 10^{-6}\f$. The function has one solution at \f$\mathbf{x} = [a, 2a]^\top\f$,
    * with \f$f(\mathbf{x}) = 0\f$. The initial guess is generated at \f$\mathbf{x} = [1, 1]^\top\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Brown : public Function<Real, 2, 3, Brown<Real>>
    {
    private:
      Real m_a{1.0e-6}; /**< Scaling value (keep it low to guarantee bad scaling). */

    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using InputVector  = typename Function<Real, 2, 3, Brown<Real>>::InputVector;
      using OutputVector = typename Function<Real, 2, 3, Brown<Real>>::OutputVector;
      using Matrix       = typename Function<Real, 2, 3, Brown<Real>>::Matrix;
      using Tensor       = typename Function<Real, 2, 3, Brown<Real>>::Tensor;

      /**
      * Class constructor for the Brown function.
      */
      Brown()
      {
        this->m_solutions.emplace_back(this->m_a, 2.0*this->m_a);
        this->m_guesses.emplace_back(1.0, 1.0);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Brown";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(const InputVector & x, OutputVector & out) const
      {
        out << x(0) - this->m_a,
               x(1) - 2.0*this->m_a,
               x(0)*x(1) - 2.0;
      }
      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(const InputVector & x, Matrix & out) const
      {
        out << 1.0, 0.0, x(1),
               0.0, 1.0, x(0);
      }

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(const InputVector & /*x*/, Tensor & out) const
      {
        out.resize(this->output_dimension());
        out[0].setZero();
        out[1].setZero();
        out[2] << 0.0, 1.0,
                  1.0, 0.0;
      }

    }; // class Brown

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_BROWN_HH
