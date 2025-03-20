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

#ifndef OPTIMIST_SCALAR_FUNCTION_COS_HH
#define OPTIMIST_SCALAR_FUNCTION_COS_HH

#include "Optimist/TestSet.hh"
#include "Optimist/Function/ScalarFunction.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the cosine function.
    *
    * Class container for the cosine function, which is defined as:
    * \f[
    * f(x) = \cos(x) \text{.}
    * \f]
    * The function has roots at \f$x = i\pi\f$, with \f$f(x) = 0\f$, and minima at \f$x = -\pi + 2i\pi\f$,
    * with \f$f(x) = -1\f$ and \f$i = 0, 1, \ldots, n\f$. The initial guesses are generated on the
    * range \f$x \in \left[-\pi, \pi\right]\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Cos : public ScalarFunction<Real, Cos<Real>>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      /**
      * Class constructor for the cosine function.
      */
      Cos()
      {
        this->m_solutions.emplace_back(-M_PI/2.0); // Zero
        this->m_solutions.emplace_back(M_PI/2.0); // Zero
        this->m_solutions.emplace_back(-M_PI); // Minimum
        this->m_solutions.emplace_back(M_PI); // Minimum
        this->m_guesses.emplace_back(-3.0/8.0*M_PI);
        this->m_guesses.emplace_back(3.0/8.0*M_PI);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Cos";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = std::cos(x);}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = -std::sin(x);}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real x, Real & out) const {out = -std::cos(x);}

    }; // class Cos

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_SCALAR_FUNCTION_COS_HH
