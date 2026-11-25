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

#ifndef OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH
#define OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH

#include "Optimist/RootFinder/OpenMethod.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Class container for the scalar Newton-Raphson method.
     *
     * \includedoc docs/markdown/RootFinder/NewtonRaphson.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class NewtonRaphson : public OpenMethod<Scalar, NewtonRaphson<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the NewtonRaphson solver.
       */
      NewtonRaphson() {}

      /**
       * Get the NewtonRaphson solver name.
       * \return The NewtonRaphson solver name.
       */
      constexpr std::string name_impl() const {
        return "NewtonRaphson";
      }

      /**
       * Newton step \f$ h = -f/f^\prime \f$.
       * \param[in] function Function (unused).
       * \param[in] x Current point (unused).
       * \param[in] f Function value.
       * \param[in] df First derivative value.
       * \param[in] ddf Second derivative value (unused).
       * \return The step.
       */
      template <typename FunctionLambda>
      Scalar step_impl([[maybe_unused]] FunctionLambda &function,
                       [[maybe_unused]] const Scalar x,
                       const Scalar f,
                       const Scalar df,
                       [[maybe_unused]] const Scalar ddf) const {
        return -f / df;
      }

    };  // class NewtonRaphson

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH
