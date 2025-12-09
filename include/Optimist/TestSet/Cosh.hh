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

#ifndef OPTIMIST_TESTSET_COSH_HH
#define OPTIMIST_TESTSET_COSH_HH

#include "Optimist/Function.hh"

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
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    requires TypeTrait<Scalar>::IsScalar
    class Cosh : public Function<Scalar, Scalar, Cosh<Scalar>>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

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
      constexpr std::string name_impl() const {return "Cosh";}

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(Scalar const x, Scalar & out) const
      {
        out = std::cosh(x);
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(Scalar const x, Scalar & out) const
      {
        out = std::sinh(x);
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(Scalar const x, Scalar & out) const
      {
        out = std::cosh(x);
        return std::isfinite(out);
      }

    }; // class Cosh

  } // namespace TestSet

} // namespace Optimist

#endif // OPTIMIST_TESTSET_COSH_HH
