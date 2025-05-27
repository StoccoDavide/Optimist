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

#ifndef OPTIMIST_TESTSET_COSH_HH
#define OPTIMIST_TESTSET_COSH_HH

#include "Optimist/TestSet.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the hyperbolic cosine function.
    *
    * Class container for the hyperbolic cosine function, which is defined as:
    * \f[
    * f(x) = \cosh(x) = \displaystyle\frac{e^x + e^{-x}}{2} \text{.}
    * \f]
    * The function has a minimum at \f$x = 0\f$, with \f$f(x) = 1\f$. The initial guesses are
    * generated on the range \f$x \in \left[-10, 10\right]\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Cosh : public Function<Real, 1, 1, Cosh<Real>>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      /**
      * Class constructor for the hyperbolic cosine function.
      */
      Cosh()
      {
        this->m_solutions.emplace_back(0.0); // Minimum
        this->m_guesses.emplace_back(-10.0);
        this->m_guesses.emplace_back(10.0);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Cosh";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = std::cosh(x);}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = std::sinh(x);}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real x, Real & out) const {out = std::cosh(x);}

    }; // class Cosh

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_COSH_HH
