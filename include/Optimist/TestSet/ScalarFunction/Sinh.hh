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

#ifndef OPTIMIST_SCALAR_FUNCTION_SINH_HH
#define OPTIMIST_SCALAR_FUNCTION_SINH_HH

#include "Optimist/TestSet.hh"
#include "Optimist/Function/ScalarFunction.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the hyperbolic sine function.
    *
    * Class container for the hyperbolic sine function, which is defined as:
    * \f[
    * f(x) = \sinh(x) = \displaystyle\frac{e^x - e^{-x}}{2} \text{.}
    * \f]
    * The function has a zero at \f$x = 0\f$, with \f$f(x) = 0\f$. The initial guesses are generated
    * on the range \f$x \in \left[-10, 10\right]\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Sinh : public ScalarFunction<Real, Sinh<Real>>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      /**
      * Class constructor for the hyperbolic sine function.
      */
      Sinh()
      {
        this->m_solutions.emplace_back(0.0); // Zero
        this->m_guesses.emplace_back(-10.0);
        this->m_guesses.emplace_back(10.0);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Sinh";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = std::sinh(x);}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = std::cosh(x);}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real x, Real & out) const {out = std::sinh(x);}

    }; // class Sinh

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_SCALAR_FUNCTION_SINH_HH
