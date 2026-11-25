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

#ifndef OPTIMIST_ROOTFINDER_CHEBYSHEV_HH
#define OPTIMIST_ROOTFINDER_CHEBYSHEV_HH

#include "Optimist/RootFinder/OpenMethod.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Class container for Chebyshev's method.
     *
     * \includedoc docs/markdown/RootFinder/Chebyshev.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Chebyshev : public OpenMethod<Scalar, Chebyshev<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{true};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Chebyshev solver.
       */
      Chebyshev() {}

      /**
       * Get the Chebyshev solver name.
       * \return The Chebyshev solver name.
       */
      constexpr std::string name_impl() const {
        return "Chebyshev";
      }

      /**
       * Chebyshev step \f$ h = -\frac{f}{f^\prime}\left(1 +
       * \frac{L}{2}\right) \f$, with \f$ L = f f^{\prime\prime} /
       * f^{\prime 2} \f$. Far from the root (\f$ |L| \geq 1 \f$) the Newton
       * step is taken.
       * \param[in] function Function (unused).
       * \param[in] x Current point (unused).
       * \param[in] f Function value.
       * \param[in] df First derivative value.
       * \param[in] ddf Second derivative value.
       * \return The step.
       */
      template <typename FunctionLambda>
      Scalar step_impl([[maybe_unused]] FunctionLambda &function,
                       [[maybe_unused]] const Scalar x,
                       const Scalar f,
                       const Scalar df,
                       [[maybe_unused]] const Scalar ddf) const {
        const Scalar newton{-f / df};
        const Scalar L{f * ddf / (df * df)};
        return std::abs(L) < 1 ? newton * (1 + L / 2) : newton;
      }

    };  // class Chebyshev

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_CHEBYSHEV_HH
