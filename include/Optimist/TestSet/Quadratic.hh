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

#ifndef OPTIMIST_TESTSET_QUADRATIC_HH
#define OPTIMIST_TESTSET_QUADRATIC_HH

#include "Optimist/TestSet.hh"

namespace Optimist
{

  namespace TestSet
  {

    /**
    * \brief Class container for the quadratic function.
    *
    * Class container for the quadratic function, which is defined as:
    * \f[
    * f(x) = ax^2 + bx + c \text{.}
    * \f]
    * The function has roots at
    * \f[
    * x = \sqrt{\displaystyle\frac{-b \pm \sqrt{b^2 - 4ac}}{2a}} \text{,}
    * \f]
    * and a minimum at
    * \f[
    * x = -\displaystyle\frac{b}{2a} \text{.}
    * \f]
    * The default coefficients are \f$a = 1\f$, \f$b = 1\f$, and \f$c = -1\f$, and the function
    * initial guess is \f$x = 0\f$.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Quadratic : public Function<Real, 1, 1, Quadratic<Real>>
    {
      Real m_a{1.0}; /**< Coefficient \f$a\f$. */
      Real m_b{1.0}; /**< Coefficient \f$b\f$. */
      Real m_c{-1.0}; /**< Coefficient \f$c\f$. */

    public:
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      /**
      * Class constructor for the quadratic function.
      */
      Quadratic()
      {
        Real delta{std::sqrt(this->m_b*this->m_b - 4.0*this->m_a*this->m_c)};
        this->m_solutions.emplace_back((-this->m_b + delta)/(2.0*this->m_a));
        this->m_solutions.emplace_back((-this->m_b - delta)/(2.0*this->m_a));
        this->m_solutions.emplace_back(this->m_b/(2.0*this->m_a));
        this->m_guesses.emplace_back(0.0);
      }

      /**
      * Get the function name.
      * \return The function name.
      */
      std::string name_impl() const {return "Quadratic";}

      /**
      * Compute the function value at the input point.
      * \param[in] x Input point.
      * \param[out] out The function value.
      */
      void evaluate_impl(Real x, Real & out) const {out = this->m_a*x*x + this->m_b*x + this->m_c;}

      /**
      * Compute the first derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The first derivative value.
      */
      void first_derivative_impl(Real x, Real & out) const {out = 2.0*this->m_a*x + this->m_b;}

      /**
      * Compute the second derivative value at the input point.
      * \param[in] x Input point.
      * \param[out] out The second derivative value.
      */
      void second_derivative_impl(Real /*x*/, Real & out) const {out = 2.0*this->m_a;}

    }; // class Quadratic

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_QUADRATIC_HH
