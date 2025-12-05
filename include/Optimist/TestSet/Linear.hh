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

#ifndef OPTIMIST_TESTSET_LINEAR_HH
#define OPTIMIST_TESTSET_LINEAR_HH

#include "Optimist/Function.hh"

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
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    requires TypeTrait<Scalar>::IsFloatingPoint
    class Linear : public Function<Scalar, Scalar, Linear<Scalar>>
    {
      Scalar m_m{1.0}; /**< Coefficient \f$ m \f$. */
      Scalar m_q{1.0}; /**< Coefficient \f$ q \f$. */

    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

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
      constexpr std::string name_impl() const {return "Linear";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Scalar x, Scalar & out) const
      {
        out = this->m_m*x + this->m_q;
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Scalar /*x*/, Scalar & out) const
      {
        out = this->m_m;
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Scalar /*x*/, Scalar & out) const
      {
        out = 0.0;
        return std::isfinite(out);
      }

    }; // class Linear

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_LINEAR_HH
