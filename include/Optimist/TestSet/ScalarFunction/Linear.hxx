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

#ifndef OPTIMIST_SCALAR_FUNCTION_LINEAR_HXX
#define OPTIMIST_SCALAR_FUNCTION_LINEAR_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the linear function.
    *
    * Class container for the linear function, which is defined as:
    * \f[
    * f(x) = mx + q \text{.}
    * \f]
    * The function has roots at
    * \f[
    * x = -\displaystyle\frac{q}{m} \text{.}
    * \f]
    * The default coefficients are \f$m = 1\f$ and \f$q = 1\f$, and the function initial guess is
    * \f$x = 0\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Linear : public ScalarFunction<Real, Linear<Real>>
    {
      Real m_m{1.0}; /**< Coefficient \f$ m \f$. */
      Real m_q{1.0}; /**< Coefficient \f$ q \f$. */

    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      /**
      * Class constructor for the linear function.
      */
      Linear()
      {
        this->m_solutions.emplace_back(-this->m_q/this->m_m);
        this->m_guesses.emplace_back(0.0);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Linear";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = this->m_m*x + this->m_q;}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real /*x*/, Real & out) const {out = this->m_m;}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real /*x*/, Real & out) const {out = 0.0;}

    }; // class Linear

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_SCALAR_FUNCTION_LINEAR_HXX
