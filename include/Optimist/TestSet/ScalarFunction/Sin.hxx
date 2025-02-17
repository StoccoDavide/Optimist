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

#ifndef OPTIMIST_SCALAR_FUNCTION_SIN_HXX
#define OPTIMIST_SCALAR_FUNCTION_SIN_HXX

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the sine function.
    *
    * Class container for the sine function, which is defined as:
    * \f[
    * f(x) = \sin(x) \text{.}
    * \f]
    * The function has roots at \f$x = \pi/2 + i\pi\f$, with \f$f(x) = 0\f$, and minima at \f$x = 3\pi/2 +
    * 2i\pi\f$, with \f$f(x) = -1\f$ and \f$i = 0, 1, \ldots, n\f$. The initial guesses are generated
    * on the range \f$x \in \left[-\pi, \pi\right]\f$.
    */
    class Sin : public ScalarFunction<Sin>
    {
    public:
      /**
      * Class constructor for the sine function.
      */
      Sin()
      {
        this->m_solutions.emplace_back(-M_PI); // Zero
        this->m_solutions.emplace_back(0.0); // Zero
        this->m_solutions.emplace_back(M_PI); // Zero
        this->m_solutions.emplace_back(-M_PI/2.0); // Minimum
        this->m_guesses.emplace_back(-3.0/4.0*M_PI);
        this->m_guesses.emplace_back(3.0/4.0*M_PI);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Sin";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = std::sin(x);}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = std::cos(x);}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real x, Real & out) const {out = -std::sin(x);}

    }; // class Sin

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_SCALAR_FUNCTION_SIN_HXX
