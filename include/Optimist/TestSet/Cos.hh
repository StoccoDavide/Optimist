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

#ifndef OPTIMIST_TESTSET_COS_HH
#define OPTIMIST_TESTSET_COS_HH

#include "Optimist/Function.hh"

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
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    requires TypeTrait<Scalar>::IsScalar
    class Cos : public Function<Scalar, Scalar, Cos<Scalar>>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the cosine function.
       */
      Cos()
      {
        this->m_solutions.emplace_back(-M_PI/2.0); // Zero
        this->m_solutions.emplace_back(0.0); // Zero
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
      constexpr std::string name_impl() const {return "Cos";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Scalar x, Scalar & out) const
      {
        out = std::cos(x);
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Scalar x, Scalar & out) const
      {
        out = -std::sin(x);
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Scalar x, Scalar & out) const
      {
        out = -std::cos(x);
        return std::isfinite(out);
      }

    }; // class Cos

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_COS_HH
