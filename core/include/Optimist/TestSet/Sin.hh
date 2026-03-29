/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_SIN_HH
#define OPTIMIST_TESTSET_SIN_HH

#include "Optimist/Function.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the sine function.
     *
     * Class container for the sine function, which is defined as:
     * \f[
     * f(x) = \sin(x) \text{.}
     * \f]
     * The function has roots at \f$x = \pi/2 + i\pi\f$, with \f$f(x) = 0\f$,
     * and minima at \f$x = 3\pi/2 + 2i\pi\f$, with \f$f(x) = -1\f$ and \f$i =
     * 0, 1,
     * \ldots, n\f$. The initial guesses are generated on the range \f$x \in
     * \left[-\pi, \pi\right]\f$.
     * \tparam Scalar Floating point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Sin : public Function<Scalar, Scalar, Sin<Scalar>> {
     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the sine function.
       */
      Sin() {
        this->m_solutions.emplace_back(-M_PI);        // Zero
        this->m_solutions.emplace_back(0.0);          // Zero
        this->m_solutions.emplace_back(M_PI);         // Zero
        this->m_solutions.emplace_back(-M_PI / 2.0);  // Minimum
        this->m_guesses.emplace_back(-0.8 * M_PI);
        this->m_guesses.emplace_back(0.8 * M_PI);
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Sin";
      }

      /**
       * Compute the function value at the input point.
       * \param[in] x Input point.
       * \param[out] out The function value.
       * \return The boolean flag for successful evaluation.
       */
      bool evaluate_impl(const Scalar x, Scalar &out) const {
        out = std::sin(x);
        return std::isfinite(out);
      }

      /**
       * Compute the first derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The first derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool first_derivative_impl(const Scalar x, Scalar &out) const {
        out = std::cos(x);
        return std::isfinite(out);
      }

      /**
       * Compute the second derivative value at the input point.
       * \param[in] x Input point.
       * \param[out] out The second derivative value.
       * \return The boolean flag for successful evaluation.
       */
      bool second_derivative_impl(const Scalar x, Scalar &out) const {
        out = -std::sin(x);
        return std::isfinite(out);
      }

    };  // class Sin

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_SIN_HH
